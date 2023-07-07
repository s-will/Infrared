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

# # RNA Design Tutorial

# ## Simple sampling of RNA sequences
# We start by importing the infrared module (and assign a shortname).

import infrared as ir

# Let's specify our first constraint network model. It is 
# _very_ simple: we define 20 variables, each with 4 possible values, which encode a sequence of nucleotides. Our first model has no dependencies.

model = ir.Model(20,4)

# From this model, we directly construct a sampler.

sampler = ir.Sampler(model)


# Once initialized with a model, the sampler can prepare the sampling from this model. In particular, it computes a tree decomposition of the constraint network. In general, it is interesting to inspect the tree decomposition of a model's constraint network, since it guides the computation. 
#
# We define a function to display information on the tree decomposition.

def show_td_info(sampler,width=600):
    td = sampler.td
    print("tree width =", td.treewidth())
    print("bags =", td.bags)
    print("edges =", td.edges)
    
    tmpfile="tmp_out.png"
    sampler.plot_td(tmpfile,'png')
    from IPython.display import Image
    return Image(filename=tmpfile,width=width)


# When we now call this function on our sampler, we will see that so far our model has a trivial tree decomposition (as expected, since we didn't add any constraints or functions to our model).

show_td_info(sampler)

# Nevertheless, we can evaluate the cluster tree.
#
# In this simple case, this will count the nucleotide sequences of length 20 (i.e. we expect to see a count of 4**20).

count = sampler.evaluate()
print("# =",int(count))

# Next, let's draw 10 samples. These samples encode uniformly drawn sequences of length 20. 

samples = [sampler.sample().values() for i in range(10)]
samples

# ... and show them more pretty.

from infrared import rna
[rna.values_to_seq(s) for s in samples]

# ## Adding complementarity constraints from an RNA secondary structure

# We define a toy RNA secondary structure (in dot-bracket notation)

structure = "((((...))))(((...)))"

# ... and parse it.

bps = rna.parse(structure)
print(bps)

# Next, we define the class of complementarity constraints; we call it
# BPComp for base pair complementarity. 
# Each such constraint is constructed for two specific positions i and j. 
# Then, the constraint is checked by testing whether the values x and y at the respective
# positions i and j correspond to complementary base pairs. 

ir.def_constraint_class( 
    'BPComp',
    lambda i,j: [i,j],
    lambda x,y: rna.values_to_seq([x,y]) 
                  in ["AU","CG","GC","GU","UA","UG"]
)    

# Btw, there is already a pre-defined constraint rna.BPComp,
# which we could have used as well.

# From the parsed structure, we generate a set of complementarity
# constraints - one for each base pair.

cons = [ BPComp( i , j ) for (i,j) in bps ]

# For illustrations, let's also print the variables of the constraints.

[c.vars() for c in cons]


# Now, we are ready to construct the new constraint model, including the
# new constraints, and construct the corresponding sampler.

# +
seqlen = len(structure) # --> number of positions / variables in the CN
model = ir.Model(seqlen,4)
model.add_constraints(cons)

sampler = ir.Sampler(model)
# -

# Let's first see whether the tree decomposition looks more interesting
# this time.

show_td_info(sampler)

# As expected, the tree decomposition reflects the dependencies due to the new constraints. 

# When we evaluate the model, we obtain the number of sequences that are
# compatbile with our secondary structure (doing simple combinatorics, we
# expect $6^7 * 6^4$).

count = sampler.evaluate()
print("# =",int(count))

# It's time to generate samples again. Before this get's boring, we quickly
# define a function to automatize this. For the fun of it, let's also draw sequence logos from the samples.
#
# **Note**: logos will be drawn only if the modules `logomaker` and `matplotlib.pyplot` are properly installed.

# +
def draw_logo(sequences):
    """Draw sequence logo for a set of sequences"""
    import logomaker as lm
    import matplotlib.pyplot as plt
        
    matrix = lm.alignment_to_matrix(sequences = sequences)
    logo = lm.Logo(matrix)
    logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    logo.ax.xaxis.set_ticks_position('none')
    plt.savefig('test.svg')
    plt.show()

def opt_draw_logo(sequences):
    """Draw sequence logo, only if required modules are available""" 
    try:
        draw_logo(sequences)
    except ModuleNotFoundError:
        pass

def spit_them_samples_out(sampler,num,structure=None):
    """Generate and show samples"""
    samples = [ sampler.sample() for i in range(num) ]
    sequences =  [ rna.ass_to_seq(s) for s in samples ]
    if structure: print(structure)
    for s in sequences: print(s)
    opt_draw_logo(sequences)


# -

# Generate samples and print them with the given structure; we can see that the base pair constraints are satisfied by all samples.

spit_them_samples_out(sampler,20,structure)

# ## Controlling the GC content
#
# To control features of our samples like the GC content, we are going to
# add suitable functions to our model. For the GC control, we define the
# function class `GCCont`, which let us count the Gs and Cs in the
# encoded sequence. For a specific position i, the function returns 1 if
# the encoded nucleotide is G or C; otherwise, 0.

ir.def_function_class(
    'GCCont',
    lambda i: [i],
    lambda x: rna.value_to_nucleotide( x ) in "GC"
)

# Btw, this function is as well pre-defined in Infrared's rna support module as
# `rna.GCCont`.
#
# GC-control is added to the sampler, simply by specifying the set of
# functions (`gc_funs`), adding them to the sampler as a function group
# `'gc'`, and setting the weight of the corresponding (automatically
# generated) feature of the same name.
#
# Note how the new functions can simply be added to the existing model with
# base pair complementarity constraints. This provides a first good example of the
# compositional construction of models in Infrared.

# +
gc_weight = 0.12 # <- try different weights: -0.1, 0, 0.1, ...

# define set of GCCont functions
gc_funs = [ GCCont( i )
              for i in range( seqlen ) ]

# add functions as group 'gc' 
model.add_functions(gc_funs, 'gc')

# set the weight of the 'gc' feature
model.set_feature_weight(gc_weight, 'gc')

# construct sampler
sampler = ir.Sampler( model )

# generate samples
spit_them_samples_out( sampler, 10, structure )
# -

# ## Controlling the BP energy
#
# Next, we are going to introduce functions to compute the energy our secondary structure 
# for the enoded sequence. To express energy in the base pair energy model, 
# we are going to define a set of functions of the pre-defined type
# `rna.BPEnergy`; one for each base pair.

bpe_funs = [ rna.BPEnergy( i, j, False ) for (i,j) in bps ] 

# These functions are added to the existing model, as the GC functions
# before. We put them into a new function group `'energy'`, which will
# immediately allow us to control them (in addition to the GC content) through the
# automatically generated feature `'energy'`.

model.add_functions(bpe_funs, 'energy')

# Let's set weights of our two features and generate samples. By setting a negative weight for energy, we aim at low energy for the structure. The GC control allows us to counter-act the induced shift towards high GC content. For the latter purpose, we set a negative weight to aim at lower GC content.

# +
model.set_feature_weight(-0.5, 'gc')
model.set_feature_weight(-1, 'energy')

sampler = ir.Sampler(model)

spit_them_samples_out(sampler, 10, structure)
# -

# ## Further constraints
# We can even extend the model by further constraints, while maintaining the ability to control the features. As example of additional hard constraints, we are going to avoid GG dinucleotides.
# Extending the model by such additional constraints, expectedly, goes through the steps of
#
# * defining the new constraint type
#
# * defining a set of constraints
#
# * adding the set to the model

# +
ir.def_constraint_class('AvoidGGConstraint',
                         lambda i: [i, i+1],
                         lambda x,y: rna.values_to_seq([x,y]) != "GG")

gg_cons = [ AvoidGGConstraint( i ) for i in range(seqlen-1) ]

model.add_constraints(gg_cons)
# -

# Finally, we construct the sampler for the extended model and generate samples.

# +
sampler = ir.Sampler(model)

spit_them_samples_out(sampler, 10, structure)
# -

# ## Targeting specific GC content and energy
#
# The previous model allowed to control the mean GC content and base pair energy of the produced model by manually tweaking the weights.
#
# This leaves several challenges if we would want to produce samples with defined GC content and energy:
#
# * we would have to filter the generated samples for a specific tolerance range around the targeted values
#
# * the features are not independent, thus changing one weight requires changing the other - this would get even harder for more than two features.
#
# Infrared supports these goals in a general way by implementing a multi-dimensional Boltzmann sampling strategy.
#
# For this purpose, first set up the model and sampler as before.

# +
model =  ir.Model(seqlen, 4)
model.add_constraints(cons)
model.add_functions(bpe_funs, 'energy')
model.add_functions(gc_funs, 'gc')

sampler = ir.Sampler( model )
# -

# Infrared now allows to set target values and tolerances for all (or selected) features and then produce targeted samples using the method `Sampler.targeted_sample()`.
#
# In the code below, we print samples with their energy and GC-content. Note how we obtain these values simply by asking the model to evaluate the features for the sample.

# +
sampler.set_target( -12, 1, 'energy' )
sampler.set_target( 10, 2, 'gc' )

samples = list()
for i in range(20):
    sample = sampler.targeted_sample()
    sequence = rna.ass_to_seq(sample)
    print("{} {:.2f} {:.2f}".format(sequence, 
                                    model.eval_feature(sample,'energy'),
                                    model.eval_feature(sample,'gc')))
    samples.append(sequence)
opt_draw_logo(samples)
# -

# ## Targeting Turner energy
#
# Finally, for real applications, it is typically much more relevant to target specific Turner energy than specific energy in the simple base pair model.
# In Infrared, we solve this by defining and adding a new feature `'Energy'` to the model that calculates Turner energy (by calling a function of the Vienna RNA package). In addition to defining this calculation, we specify as well that the new feature should control the function group `energy`. In this way, Infrared will use base pair energy as proxy of Turner energy.
#
# **Note:** this requires the Vienna RNA package with working Python bindings (currently, this fails in Windows even after installing the package from binaries)

# +
import RNA

model.add_feature( 'Energy', # feature name
                   'energy', # controlled group(s)
                   #
                   # function to evaluate the feature for a sample;
                   # NOTE how we have to bind i
                   lambda sample, structure=structure:
                      RNA.energy_of_struct( rna.ass_to_seq( sample ),
                                            structure )
                 )
# -

# As before, we set targets and generate samples. Of course, this time, we set a target range for the new feature of Turner energy.

# +
sampler = ir.Sampler(model)

sampler.set_target( -5, 1, 'Energy' )
sampler.set_target( 10, 2, 'gc' )

samples = list()
for i in range(20):
    sample = sampler.targeted_sample()
    print("{} {:5.2f} {:5.2f} {:5.2f}".format(rna.ass_to_seq(sample), 
                                    model.eval_feature(sample,'energy'),
                                    model.eval_feature(sample,'Energy'),
                                    model.eval_feature(sample,'gc')))
    samples.append(rna.ass_to_seq(sample))

opt_draw_logo(samples)
# -

# ## Add IUPAC constraints
#
# In many design applications, we have prior knowledge on the sequences, which can be encoded as IUPAC string. 
# Below, we generate constraints from a IUPAC string and add them to the Infrared model for targeting Turner energy and GC content.

# +
sequence = "RSSSUWWSSNNSNNNNMNYR"

for i,x in enumerate(sequence):
    model.add_constraints(ir.ValueIn(i,rna.iupacvalues(x)))

sampler = ir.Sampler(model)

sampler.set_target( -5, 1, 'Energy' )
sampler.set_target( 10, 2, 'gc' )

samples = list()
for i in range(20):
    sample = rna.ass_to_seq(sampler.targeted_sample())
    print(sample)
    samples.append(sample)
    
opt_draw_logo(samples)
# -

# **Remark:** At this point, we have implemented the full functionality of the RNA design approach IncaRNAtion (Reinharz et al., Bioinformatics 2013) in Infrared. In Infrared, we can easily extend the model further by including further constraints and going on to multi-target design. Infrared makes the latter, which extends functionality to the tool RNARedPrint (Hammer et al; BMC Bioinf 2021), look surprisingly simple. (We demonstrate this in a separate accompanying notebook, as well as in the Infrared-based application RNARedPrint 2).
