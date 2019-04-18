#!/usr/bin/env python3
# This file is a plug-in of the InfraRed to realise the sampling of module sequences using 
# BayesPairing

import pickle
import networkx as nx
import infrared as ir

# TODO: Modify the function to return the proper value
class LogProb(ir.Function):
    def __init__(self, lst, p):
        self.p = p
        self.index = lst

    def vars(self):
        return self.index
        
    def __call__(self, seq):
        prob = -10
        return prob

class ModuleConstraintNetwork(ir.ConstraintNetwork):
    def __init__(self, module):
        self.dependencies = list(map(lambda t: list(module.predecessors(t))+[t], module.nodes())) 

        self.functions = list(map(lambda t: (t, [LogProb([],1)]),self.dependencies))
        self.constraints = []
        self.seqlen = len(module.nodes())

class ModuleTreeDecomposition(ir.TreeDecomposition):
    def __init__(self, cn, *, method=0):
        super().__init__(cn.seqlen,cn.dependencies, method=method)
        self.domains = 4
        self.cn = cn

    def get_bag_assignments(self):
        print("Start Assign Functions")
        bagfunctions = self.assign_to_bags(self.cn.functions)
        bagconstraints = self.assign_to_bags(self.cn.constraints)
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

def main(f):
    module = pickle.load(open(f, 'rb'))
    new_module = nx.convert_node_labels_to_integers(module)
    sampler = ModuleSampler(new_module)
    sampler_generator = sampler.samples()
    return sampler_generattackr

if __name__ == "__main__":
    main("bayes_net_hairpin.cPickle")


