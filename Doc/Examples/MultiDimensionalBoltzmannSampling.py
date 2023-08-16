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

# # Multi-dimensional Boltzmann sampling
#
# We demonstrate sampling of sequences targeting specific dinucleotide frequences and additional requirements; 
# in this example, compatibility with an RNA target structure.
#
# The targeting of dinucleotide frequencies is performed using the multi-dimensional Boltzmann sampling functionality of Infrared's Sampler (aka MultiDimensionalBoltzmannSampler) class.
#

# ## Import modules + define timing function

# +
import infrared as ir
import infrared.rna as rna
from collections import Counter
import random

import time
def report_time(f,desc=""):
    start=time.time()
    res=f()
    desc = "" if desc=="" else f' ({desc})'
    print(f'Time{desc}: {time.time()-start:0.2f}s')
    return res


# -

# ## Preliminaries

alphabet = "ACGU"
dimers = [ x+y for x in alphabet for y in alphabet ]

# ## Network function for dinucleotide counting
#
# For use in the model, define counting of dinucleotide frequencies. The generated class will be constructed by the position $i$ (to count the dimer at positions $i$, $i+1$) and the respective dimer as string ("AA", "AC", ...).
#
# Note how the dimer argument is passed from the first lambda function (`init'), which specified the construction and returns the dependency list [$i,i+1$] , to the second one (`value'), which evaluates the function given the values of variables $i$ and $i+1$.

ir.def_function_class('DinuclFreq',
    lambda i,dimer: [i,i+1],
    lambda x,y,dimer: rna.values_to_seq([x,y])==dimer
)

# ## Input instance(s)

# +
#>CP000576.1:c381532-381475 Prochlorococcus marinus str. MIT 9301, complete genome
#UUUCGUUCACCCUCAUUUGAGGGCGCAGUUCGAGUCAUACCAUGGAACGGGGGAUGGC
seq0    = "UUUCGUUCACCCUCAUUUGAGGGCGCAGUUCGAGUCAUACCAUGGAACGGGGGAUGGC"
struct0 = "...[[[[[.(((((....)))))..........(((((.((...]]]]]))..)))))"

## Yann's example
seq = "CACUGUCGACUCAGUCAGUGAGUCGCGACUGACUGCAUCGCGACUACGUCAGUCGUACGAUCUAUUGUGCGAUAUCGCGCGAUAUAGCUAUGCG"
struct = None

# SAM (RF00162)
#>AP006840.1:c2688754-2688649 Symbiobacterium thermophilum IAM 14863 DNA, complete genome
#GGUUCAUCGAGAGUGGCGGAGGGACUGGCCCCAUGAUGCCACGGCAACCUCUCCCGCGGGGAGAACGGUGCCAAAUCCAGCGGACACUCGGUCCGAGAGAUGAAGC
seq    = "GGUUCAUCGAGAGUGGCGGAGGGACUGGCCCCAUGAUGCCACGGCAACCUCUCCCGCGGGGAGAACGGUGCCAAAUCCAGCGGACACUCGGUCCGAGAGAUGAAGC"
struct = "(.((((((....((((((..(((.[[[[.)))....))))))((((.((((((((...))))))..))))))....]]]](((((.....)))))...)))))).)"
print("Sequence length",len(seq))
# -

# Obtain target frequencies from the sequence; print the counts

tgtfreqs = Counter(seq[i]+seq[i+1] for i in range(len(seq)-1))
print(tgtfreqs)

# ## Feature network model for the design task

# ## Generate model
#
# The model consists of
#
# * $n$ (sequence length) many variables, each with domains 0..3 (enconding the nucleotides);
#
# * constraints to guarantee compatibility with the given structure;
#
# * functions counting the different dinucleotides; defining one feature per dimer.

# +
n = len(seq)
model = ir.Model(n,4)

# complementarity constraints
if struct is not None:
    model.add_constraints(rna.BPComp(*bp) for bp in rna.parse(struct))

# dinucleotide frequency functions
for dimer in dimers:
    model.add_functions([DinuclFreq(i,dimer) for i in range(n-1)],dimer)
# -

# ## Draw samples (without targeting)
#
# This is mainly done, to check run-times for generating a number of samples (as is done by the multi-dimensional Boltzmann sampling strategy in every iteration). As well report the treewidth of the network (as determined by the sampler).

nsamples=100
sampler = report_time(lambda:ir.Sampler(model,lazy=False),'Construction')
print("Treewidth:",sampler.treewidth())
_ = report_time(lambda:[rna.ass_to_seq(sampler.sample()) for _ in range(n)],'Sampling')

# ## Draw targeted samples

# Generate sampler, set parameters for the multi-dimensional Boltzmann sampling strategy and set the targets.

# +
sampler = ir.Sampler(model,lazy=True)
sampler.samples_per_round = 200
sampler.tweak_factor = 0.05

tolerance = 2

for dimer in dimers:
    sampler.set_target(tgtfreqs[dimer],tolerance=tolerance,featureid=dimer)


# +
## collect some statistics
def callback(total,accepted,fstats):
    rmsd = sampler.rmsd(fstats.means(),sampler.model.features)
    the_statistics.append({'samples':total, 'rmsd':rmsd, 'accepted':accepted})
    #print(fstats.report(),rmsd)
    
sampler.callback = callback
# -

# Actually, draw the targeted samples (and stop run time.)

# +
sampler.verbose = False

the_statistics = list()
the_samples = list() ## store the produced samples

print_samples = None

stats = ir.FeatureStatistics()
def draw_samples(n):
    for i in range(n):
        sample = sampler.targeted_sample()
        feature_values = {dimer:model.eval_feature(sample,dimer) for dimer in dimers}
        stats.record_features(model.features,feature_values)

        if print_samples:
            line = rna.ass_to_seq(sample)
            for dimer,freq in feature_values.items(): line += f" {dimer}:{freq}"
            if type(print_samples)!=int or (i+1)%print_samples==0: 
                print(i+1, line)
        the_samples.append(sample)
        
report_time(lambda:draw_samples(1000))
# -

# Finally, report some statistics and learned weights

print("Targets",[f"{dimer}:{tgtfreqs[dimer]}" for dimer in dimers])
print("Stats  ",stats.report())
print("Weights",{k:f'{f.weight:.3f}' for k,f in sampler.model.features.items()})    

# ## Plots sampling effectivity 

# +
import matplotlib.pyplot as plt
import seaborn as sns

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(2.5,5))

sns.scatterplot(x=[row['samples'] for row in the_statistics], y=[100*row['accepted']/row['samples'] for row in the_statistics],ax=ax1)
ax1.set_ylabel("Accepted percentage")
#sns.scatterplot(x=[row['samples'] for row in the_statistics], y=[row['accepted'] for row in the_statistics],ax=ax1)
#ax1.set_ylabel("Accepted samples")
sns.scatterplot(x=[row['samples'] for row in the_statistics], y=[row['rmsd'] for row in the_statistics],ax=ax2)
ax2.set_xlabel("Total samples")
ax2.set_ylabel("RMSD")
fig.tight_layout()
plt.savefig('mdbs_run.pdf')
plt.show()
# -

# ## Distribution plots

import math
def distribution_heatmaps(model,n,feature_list,fig,*,limits,labels=None,targets=None,ax=None,cmap="Blues"):
    sampler = ir.Sampler(model)
    samples = [sampler.sample() for _ in range(n)]
    
    def eval_features(sample):
        return {f:model.eval_feature(sample,f) for f in feature_list}
    
    features = [eval_features(sample) for sample in samples]

    if labels is None:
        labels = feature_list
    
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
            #sharex=True, sharey=True)
            
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
        pad=0.4,
        location="left",
        boundaries=[0.2+i/10 for i in range(9)],
        values=[i/8 for i in range(1,9)])
    
    return ax


# +
sel_dimers = dimers #["CC","GU", "GG", "UU"]

limits=[(-0.5,18)]*len(sel_dimers)

n = 500

k=len(sel_dimers)

#fig = plt.figure(figsize=(10,5))
fig = plt.figure(figsize=(20,16))

ax = distribution_heatmaps(model,n,sel_dimers,fig,limits=limits,cmap="Blues")
distribution_heatmaps(sampler.model,n,sel_dimers,fig,limits=limits,
                      targets=tgtfreqs,
                      ax=ax,cmap="Reds")
fig.tight_layout()
plt.savefig('mdbs_heatmaps.pdf')
plt.show()
# -

