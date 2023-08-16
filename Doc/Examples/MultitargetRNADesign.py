# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Multi-target RNA Design
#
# Implements RNARedprint-like design of RNA sequences with very specific properties towards multiple targets (several target secondary structures and GC content.
#
#
# **Software requirements**
#
# Software can be installed using ```conda``` and ```pip``` on MacOSX or Linux.
#
# Apart from, obviously, ```Infrared``` (*install* by ```conda install -c conda-forge infrared```), this notebook requires
# * Vienna RNA package,
#   *install* by ```conda install -c conda-forge -c bioconda viennarna```
#
# The following software/Python modules are optionally used to generate graphical output:
# * seaborn
#   e.g. *install* via ```conda install seaborn```
#
# * VARNA and inkscape to draw target structures: 
#   *install* VARNA (and its Python api) by 
#      * ```pip install varnaapi```
#      * then, download VARNA (VARNAv3-92.jar or newer, https://varna.lri.fr/bin/VARNAv3-92.jar)
#        and put it into your working directory.
#      * ```conda install inkscape```
#      
# * logomaker to draw sequence logos, *install* by ```conda install logomaker```

# +
import sys
import subprocess
from IPython.display import Image, display, SVG
import tempfile
from collections import defaultdict

import infrared as ir
from infrared import rna
import RNA

try:
    import seaborn as sns
    import matplotlib.pyplot as plt
except:
    print("Get seaborn to draw sample histograms.")

try:
    import logomaker
except:
    print("Get logomaker to draw sequence logos.")

try:
    import varnaapi
except:
    print("Get VARNA to plot target structures.")

import time
def report_time(f,desc=""):
    start=time.time()
    res=f()
    desc = "" if desc=="" else f' ({desc})'
    print(f'Time{desc}: {time.time()-start:0.2f}s')
    return res
    
def has_module(m): return m in sys.modules


# -

# ## Target structures
#
# RNA target structures are specified in the form of dot bracket strings (as e.g. used by the Vienna RNA package).
# Dot-bracket strings are linear representations of RNA secondary structure, mainly used to annotate RNA sequences. They contain one character per sequence position. Each pair of balanced opening and closing brackets specifies one base pair of the secondary structure.
#
# We specify some examples and plot them using VARNA (if it is available).

# +
## our target RNA secondary structure
targets = [
#    012345678901234567890
    "((((...)))).(((...)))",
    "((((((......)))...)))",
    "......(((...)))......"
]

def plotRNA(structure, sequence=None, filename=None, show=True):
    n=len(structure)
    if sequence is None:
        sequence = ' '*n
    
    if filename is not None and not filename.endswith(".svg"):
       filename+=".svg" 

    f = filename
    if filename is None:
        f = tempfile.NamedTemporaryFile(suffix=".svg")
        filename = f.name
    
    v = varnaapi.Structure(structure=structure,sequence=sequence)
    v.set_algorithm('radiate')
    v.savefig(filename)
    subprocess.run(["inkscape", "-D", filename, "--export-overwrite"], capture_output=True)
    if show:
        display(SVG(filename))
    
    return f
    
for i, target in enumerate(targets):
    print(f"Target {i}: {target}")
    if has_module("varnaapi"):
        plotRNA(target)


# -

# ## Construct the Infrared model

def construct_design_model(targets):
    seqlen = len(targets[0])
    ######
    # construct the constraint model
    model = ir.Model()

    # one variable X_i per position i;
    # the value of X_i encodes the nucleotide at position i   
    model.add_variables( seqlen, 4 )


    for i,target in enumerate(targets):
        bps = rna.parse(target)

        model.add_constraints( rna.BPComp( i, j ) for ( i, j ) in bps )

        model.add_functions( [ rna.BPEnergy( i, j, False ) 
                               for ( i, j ) in bps ], group = f'bpenergy{i}' )

        model.add_feature( f'E{i}', # feature name
                           f'bpenergy{i}', # controlled group(s)
                           #
                           # function to evaluate the feature for a sample;
                           # NOTE how we have to bind i
                           lambda sample, i=i: RNA.energy_of_struct( rna.values_to_seq( sample.values() ),
                                                  targets[i] )
                         )

    model.add_functions( [ rna.GCCont( i = i ) for i in range(seqlen) ], group = 'gc' )


    model.write_graph('dependency_graph.dot', True)
    ir.dotfile_to_pdf('dependency_graph.dot')

    # the model generates automatic features 'bpenergyI', 'gc' from the function groups;
    # as well as total feature combining all function groups;
    # however, we want to diretly control Turner energy (instead of base pair energy).
    # For this purpose, add additional features 'EI'

    return model


# ## Some helper functions for showing the outcome of sampling

# +
def draw_logo(sequences, filename=None):
    """Draw sequence logo for a set of sequences"""
    import logomaker as lm
    import matplotlib.pyplot as plt
        
    matrix = lm.alignment_to_matrix(sequences = sequences)
    logo = lm.Logo(matrix)
    logo.style_xticks(rotation=90, fmt='%d', anchor=0, spacing=5)
    logo.ax.xaxis.set_ticks_position('none')
    if filename is not None: plt.savefig(filename)
    plt.show()

def print_sample(design):
    seq = rna.values_to_seq( design.values() )
    print(seq)

def show_designs(designs, maxdesigns=10, plotlogo=False, plotdist=False, returnresults=False, filename=None):
    """
    Produce report of sequence designs
    Args:
        designs    list of assignments of the design model
        maxdesigns show first maxdesigns designs
        plotlogo   if true, plot logo
    """
    
    sequences = []
    
    statistics = defaultdict(list)
    
    for i,design in enumerate(designs):
        seq = rna.values_to_seq( design.values() )
        sequences.append(seq)
        
        
        stats = { f'E{j}': model.eval_feature(design, f'E{j}') for j in range(len(targets))}
        stats['GC'] = round(model.eval_feature(design, 'gc')*100/len(seq))

        for k in stats:
            statistics[k].append(stats[k])
        
        if i<maxdesigns:
            print(seq,end="")
            for k in stats:
                print(f" {k}={stats[k]:.2f}",end="")
            print()
        
    if maxdesigns<len(designs):
        print(f"... skip {len(designs)-maxdesigns} designs")
    
    if plotlogo:
        print()
        print(f"Sequence logo from {len(designs)} sequences:")
        draw_logo(sequences, 
                  filename=filename+"_logo.svg" if filename is not None else None)

    
    if plotdist:
        fig,axs = plt.subplots(1, 2, figsize=(10, 2.5), width_ratios=[1,2])

        # plot only gc in first histogram
        sns.histplot(statistics['GC'],ax=axs[0], stat='probability', discrete=True, color=[0.8,0.2,0,0.6])
        
        del statistics['GC'] # take out gc for second histogram
        
        # make Ei one based for the plot
        statistics = { f'E{i+1}': statistics[f'E{i}'] for i in range(len(statistics)) }

        sns.histplot(statistics, ax=axs[1], stat='probability',
                     palette={'E1':[0.0,0.8,0.2,0.5],
                              'E2':[0.0,0.3,0.9,0.5],
                              'E3':[0.8,0.4,0.1,0.5]}
                    )
        axs[1].set_ylabel("")
        
        if filename is not None:
            fig.savefig(filename+"_hist.svg")
        fig.show()
        
    if returnresults:
        return sequences, statistics
# -

# ## Sampling from the design model

# +
model = construct_design_model(targets)

model.set_feature_weight( -2, 'bpenergy0' )
model.set_feature_weight( -2, 'bpenergy1' )
model.set_feature_weight( -2, 'bpenergy2' )
model.set_feature_weight( -1, 'gc' )

## Sampling at specific weights

sampler = ir.BoltzmannSampler( model )

print( "Tree width:", sampler.treewidth() )
print()

#opionally, write plot of tree decomposition to file
#tdfile = "treedecomp.pdf"
#sampler.plot_td(tdfile)

samples = report_time(lambda:[sampler.sample() for _ in range(1000)])

show_designs(samples,
    plotlogo=has_module("logomaker"),
    plotdist=has_module("matplotlib"))
# -

# ## Targeting features by Multi-dimensional Boltzmann sampling

# +
print("###########################################")    
## MDBS

model.set_feature_weight( 0, 'E0' )
model.set_feature_weight( 0, 'E1' )
#model.set_feature_weight( 0, 'E2' )
model.set_feature_weight( 0, 'gc' )

######
# create sampler
sampler = ir.Sampler( model )

######
# set targets

# control number of gc's; we target 70% +/- 15% GC-content
seqlen = len(targets[0])
sampler.set_target( 0.85 * seqlen, 0.02 * seqlen, 'gc' )

# control Turner energy, target -2 +/- 1 kcal/mol
sampler.set_target( -2, 0.2, 'E0' )

# control Turner energy, target -2 +/- 1 kcal/mol
sampler.set_target( -3, 0.2, 'E1' )

# control Turner energy, target -2 +/- 1 kcal/mol
#sampler.set_target( -1.5, 0.2, 'E2' )

######
# and print samples
for i in range(10):
    sample = sampler.targeted_sample()
    print_sample(sample)
# -

# # APPENDIX

# +
# RNAfold/3str/f3.100.0.inp

benchmark_targets = [
'((((.((....)).)))).((.(((.((((.....(((..((((((.((..(((.(.....).)))..)).)).))))..)))..)))).))).))....',
'..(((((.....(((.(((((((.....))))..))).))).....)))))..((((((((((...))).)....))))))...((((((....))))))',
'......(((((.....(((...(((.((.((.(((....((......))...))).)).)))))..))).............))))).((((...)))).'
]



# +
model = construct_design_model(benchmark_targets)

for weight in [-5,0]:
    for i in range(3):
        model.set_feature_weight( weight, f'bpenergy{i}' )
    model.set_feature_weight( weight, 'gc' )

    ## Sampling at specific weights

    sampler = ir.BoltzmannSampler( model )

    samples = report_time(lambda:
        [sampler.sample() for _ in range(5000)])
    
    show_designs(samples,
        plotlogo=has_module("logomaker"),
        plotdist=has_module("matplotlib"),
        filename=f"f3_{weight}")


# -

def model_objective(model, assignment):
    return sum(f.weight*model.eval_feature(assignment,fname) for fname, f in model.features.items() if f.weight!=0)


# +
model = construct_design_model(benchmark_targets)
model.set_feature_weight( -3, 'E0' )
model.set_feature_weight( -4, 'E1' )
model.set_feature_weight( -5, 'E2' )
model.set_feature_weight( -6, 'gc' )

sampler = ir.Sampler(model)
# draw samples
for i in range(10):
    sample = sampler.sample()
    print_sample(sample)
    print([(f.weight, model.eval_feature(sample,fname)) for fname, f in model.features.items() if f.weight!=0])
    print(model_objective(model,sample))


# -

# ## MDBS targeting distribution plots

# +
## collect some statistics
def callback(total,accepted,fstats):
    rmsd = sampler.rmsd(fstats.means(),sampler.model.features)
    the_statistics.append({'samples':total, 'rmsd':rmsd, 'accepted':accepted})
    #print(fstats.report(),total,accepted,rmsd)
    
the_statistics = list()


# +
model = construct_design_model(benchmark_targets)

model.set_feature_weight( 0, 'E0' )
model.set_feature_weight( 0, 'E1' )
model.set_feature_weight( 0, 'E2' )
model.set_feature_weight( 0, 'gc' )

######
# create sampler
sampler = ir.Sampler( model )

sampler.samples_per_round = 500
sampler.tweak_factor = 0.1
sampler.cooling = 2**(1/16)
sampler.verbose = True

sampler.callback = callback


######
# set targets

# control number of gc's
seqlen = len(benchmark_targets[0])
sampler.set_target( 0.85 * seqlen, 0.05 * seqlen, 'gc' )

# control Turner energy
sampler.set_target( -40, 0.5, 'E0' )

# control Turner energy
sampler.set_target( -40, 0.5, 'E1' )

# control Turner energy
sampler.set_target( -30, 0.5, 'E2' )

######
# and generate targetd samples
the_samples = report_time(
    lambda: [sampler.targeted_sample() for i in range(100)])

# +
import matplotlib.pyplot as plt
import seaborn as sns

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(2.5,5))

sns.scatterplot(x=[row['samples'] for row in the_statistics], y=[100*row['accepted']/row['samples'] for row in the_statistics],ax=ax1)
ax1.set_ylabel("Accepted percentage")
sns.scatterplot(x=[row['samples'] for row in the_statistics], y=[row['rmsd'] for row in the_statistics],ax=ax2)
ax2.set_xlabel("Total samples")
ax2.set_ylabel("RMSD")
fig.tight_layout()
plt.savefig('mdbs_design_run.pdf')
plt.show()

#print(the_statistics[-1])
# -

import math
def distribution_heatmaps(model,n,feature_list,fig,*,limits,labels,targets=None,ax=None,cmap="Blues"):
    sampler = ir.Sampler(model)
    samples = [sampler.sample() for _ in range(n)]
    
    def eval_features(sample):
        return {f:model.eval_feature(sample,f) for f in feature_list}
    
    features = [eval_features(sample) for sample in samples]

    k = len(feature_list)
    
    the_plots = [(i,j) for i in range(1,k) for j in range(0,i)]
    nplots = len(the_plots)
    dimx = int(math.sqrt(nplots))
    dimy = nplots // dimx
    
    the_cells = [(i,j) for i in range(0,dimx) for j in range(0,dimy)]
    
    if ax is None:
        ax = fig.subplots(dimx,dimy+1,
                          squeeze=False,
                          width_ratios=[1]*dimy+[1]
                         )

    for x in range(dimx):
        for y in range(dimy+1):
            if (x,y) not in the_cells[:nplots]:
                ax[x][y].axis("off")

    for idx,(i,j) in enumerate(the_plots):

        fi = feature_list[i]
        fj = feature_list[j]
        
        cax = ax[the_cells[idx][0]][the_cells[idx][1]]
        
        sns.kdeplot(x=[f[fj] for f in features],y=[f[fi] for f in features], ax = cax,
                cmap=cmap, levels=8, thresh=0.2, cbar = False, fill=True)
        if targets:
            cax.axhline(y = targets[fi], color = 'red', linestyle = 'dashed')
            cax.axvline(x = targets[fj], color = 'red', linestyle = 'dashed')

        cax.set_ylabel(labels[i])
        cax.set_xlabel(labels[j])

        cax.set_ylim(limits[i])
        cax.set_xlim(limits[j])

    plt.colorbar(plt.cm.ScalarMappable(cmap=cmap),ax=ax[0][dimy],
        pad=0.2,
        location="left",
        boundaries=[0.2+i/10 for i in range(9)],
        values=[i/8 for i in range(1,9)])
    
    return ax


# +
sel_features = ["E0","E1","E2","gc"]
labels = ["E1","E2","E3","GC content"]
limits = [(-50,15),(-50,15),(-50,15),(30, seqlen)]


n = 1000

k=len(sel_features)

targets = {k:f.target for k,f in sampler.model.features.items() if k in sel_features}


fig = plt.figure(figsize=(10,5))
ax = distribution_heatmaps(model,n,sel_features,fig,
                           limits=limits,labels=labels,cmap="Blues")
distribution_heatmaps(sampler.model,n,sel_features,fig,
                      limits=limits,labels=labels,
                      targets=targets,
                      ax=ax,cmap="Reds")
fig.tight_layout()
plt.savefig('mdbs_design_heatmaps.pdf')
plt.show()
# -

