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
sampler.samples_per_round = 100
sampler.tweak_factor = 0.025
#sampler.verbose = True

tolerance = 2

for dimer in dimers:
    sampler.set_target(tgtfreqs[dimer],tolerance=tolerance,featureid=dimer)
# -

# Actually, draw the targeted samples (and stop run time.)

# +
stats = ir.FeatureStatistics()

def draw_samples(n):
    for i in range(n):
        sample = sampler.targeted_sample()
        feature_values = {dimer:model.eval_feature(sample,dimer) for dimer in dimers}
        stats.record_features(model.features,feature_values)

        line = rna.ass_to_seq(sample)
        for dimer,freq in feature_values.items(): line += f" {dimer}:{freq}"
        print(i+1, line)

report_time(lambda:draw_samples(20))
# -

# Finally, report some statistics and learned weights

print("Targets",[f"{dimer}:{tgtfreqs[dimer]}" for dimer in dimers])
print("Stats  ",stats.report())
print("Weights",{k:f'{f.weight:.3f}' for k,f in sampler.model.features.items()})    


