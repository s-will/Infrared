#!/usr/bin/env python3
# This file is a plug-in of the InfraRed to realise the sampling of module sequences using 
# BayesPairing

import random
import pickle
import math
import networkx as nx
import pandas as pd
import argparse
from pgmpy.models import BayesianModel
import infrared as ir

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
        return p
        
class ModuleConstraintNetwork(ir.ConstraintNetwork):
    def __init__(self, seqlen, structure, modules, features):
        super().__init__()
        self.seqlen = seqlen
        self.features = features
        self.functions = list()
        module_dependencies = list()
        module_edges = set()
        
        for module in modules:
            dependencies, probs = relabel_rna_module(module)
            module_edges.update(set(module.edges()))
            self.functions.extend([(t[0], [LogProb(*t)]) for t in zip(dependencies, probs)])
            module_dependencies.extend(dependencies)
        structureset = set(structure)
        bp_dependencies = [ [i,j] for (i,j) in (structureset - module_edges)]
        self.dependencies = module_dependencies + bp_dependencies    

        self.functions.extend( [ ([i,j], [ir.BPEnergy( i,j, not (i-1,j+1) in structureset, 1)])
            for [i,j] in bp_dependencies ] )

        self.constraints = [([i,j], [ir.ComplConstraint(i,j)]) for [i,j] in bp_dependencies]
        
class ModuleTreeDecomposition(ir.TreeDecomposition):
    def __init__(self, cn, *, method=0):
        super().__init__(cn.seqlen, cn.dependencies, method=method)
        self.domains = 4
        self.cn = cn

    def get_bag_assignments(self):
        bagfunctions = self.assign_to_bags(self.cn.functions)
        bagconstraints = self.assign_to_all_bags(self.cn.constraints)
        return bagconstraints, bagfunctions

## @brief Sampler for module
class ModuleSampler(ir.BoltzmannSampler):
    def __init__(self, seqlen, structure, modules, method = 0, features = []):
        super().__init__(features)
        self.seqlen = seqlen
        self.modules = modules
        self.method = method
        self.structure = parseRNAStructureBps(structure)
        self.setup_engine()

    def gen_constraint_network(self, features):
        return ModuleConstraintNetwork(self.seqlen, self.structure, self.modules, features)
    def gen_tree_decomposition(self, cn):
        return ModuleTreeDecomposition(cn, method = self.method)

    def sample(self):
        return ir.values_to_sequence(super().sample().values())

def read_inp(f):
    structure = next(f).rstrip('\n')
    modules = list()
    for line in f:
        lst = line.rstrip('\n').split()
        m = pickle.load(open(lst[-1], 'rb'))
        nodes = list(m.nodes())
        m.map = {nodes[i]: int(lst[i]) for i in range(len(nodes))}
        modules.append(m)
    return structure, modules

def parseRNAStructureBps(structure, *, opening = "([{<", closing = ")]}>"):
    stack = { op:list() for op in opening }
    lst = [-1]*len(structure)

    for i,c in enumerate(structure):
        for (op,cl) in zip(opening,closing):
            if c==op:
                stack[op].append(i)
            elif c==cl:
                j = stack[op].pop()
                lst[i] = j
                lst[j] = i
    bps = list()
    for i,j in enumerate(lst):
        if i != -1 and i<j:
            bps.append( (i,j) )
    return bps

# #############################################################
# Main

def main(args):
    # Init seed
    if args.seed == None:
        args.seed = random.randint(0, 2**31)
    ir.seed(args.seed)
    print(args.seed)

    ir.set_bpenergy_table(ir.params_bp)

    # Read input file
    with open(args.infile) as infh:
        structure, modules = read_inp(infh)

    # Variables number
    seqlen = len(structure)
    
    sampler = ModuleSampler(seqlen, structure, modules, method = args.method)
    sampler_generator = sampler.samples()

    for count in range(args.number):
        seq = next(sampler_generator)
        print(seq)

## @brief Relabel dependencies and cpds for a given module with index map
def relabel_rna_module(module):
    cpds = module.get_cpds()
    variables = [t.variables for t in cpds]
    dependencies = [[module.map[a] for a in t] for t in variables]
    probs = [t.get_values() for t in cpds]
    return dependencies, probs

if __name__ == "__main__":
    ## argument parser
    parser = argparse.ArgumentParser(description="Boltzmann sampling for RNA design with BayesPairing modules")
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--method', type=int, choices=range(6), default=0,
        help="Method for tree decomposition")
    parser.add_argument('-n', '--number', type=int, default=10, help="Number of samples")
    parser.add_argument('--seed', type=int, default=None, help="Seed for infrared")
    args = parser.parse_args()
    main(args)
