#!/usr/bin/env python3
# This file is a plug-in of the InfraRed to realise the sampling of module sequences using 
# BayesPairing

import pickle
import math
import networkx as nx
import infrared as ir

# TODO: Modify the function to return the proper value
class LogProb(ir.Function):
    def __init__(self, indices, probs):
        super().__init__(indices)
        self.indices = indices
        self.len = len(indices)
        self.probs = probs

    # Overwrite
    def __call__(self, a):
        a = a.values()
        sub_a = [a[i] for i in self.indices]
        col = 0
        for ind, v in enumerate(sub_a[-1:0:-1]):
            col += v*(4**ind)
        p = self.probs[sub_a[0],col]
        print(p)
        return -math.log(p)
        
class ModuleConstraintNetwork(ir.ConstraintNetwork):
    def __init__(self, module):
        super().__init__()
        self.map_index = module.map
        dependencies, probs = relabel_rna_module(module)
        def to_index(lst):
            return list(map(lambda t: self.map_index[t], lst))
        self.dependencies = dependencies
        self.functions = [(t[0], [LogProb(*t)]) for t in zip(dependencies, probs)]
        self.constraints = []
        self.seqlen = len(module.nodes())

class ModuleTreeDecomposition(ir.TreeDecomposition):
    def __init__(self, cn, *, method=0):
        super().__init__(cn.seqlen, cn.dependencies, method=method)
        self.domains = 4
        self.cn = cn

    def get_bag_assignments(self):
        print("Start Assign Functions")
        bagfunctions = self.assign_to_bags(self.cn.functions)
        bagconstraints = self.assign_to_all_bags(self.cn.constraints)
        return bagconstraints, bagfunctions

## @brief Sampler for module
class ModuleSampler(ir.BoltzmannSampler):
    def __init__(self, module, method = 0, features = []):
        super().__init__(features)
        self.seqlen = module.size()
        self.module = module
        self.method = method

        self.setup_engine()

    def gen_constraint_network(self, features):
        return ModuleConstraintNetwork(self.module)

    def gen_tree_decomposition(self, cn):
        return ModuleTreeDecomposition(cn, method = self.method)

    def sample(self):
        return ir.values_to_sequence(super().sample().values())

def main(map_index, f):
    module = pickle.load(open(f, 'rb'))
    module.map = map_index
    sampler = ModuleSampler(module)
    sampler_generator = sampler.samples()
    count  = 0
    for seq in sampler_generator:
        print(seq)
        count += 1
        if count == 10:
            break

def relabel_rna_module(module):
    cpds = module.get_cpds()
    variables = [t.variables for t in cpds]
    dependencies = [[module.map[a] for a in t] for t in variables]
    probs = [t.get_values() for t in cpds]

    return dependencies, probs
    
if __name__ == "__main__":
    main({1797:0, 1798:1, 1799:2, 1800:3, 1801:4, 1802:5, 1803:6}, "test_module")
