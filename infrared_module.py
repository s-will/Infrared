#!/usr/bin/env python3

import itertools
import infrared as ir

import treedecomp
import rnastuff as rna

def val2nucl(x):
    return "ACGU"[x]
def values2seq(xs):
    return "".join(map(val2nucl,xs))

def seed(x):
    ir.seed(x)

def set_bpenergy_table(tab):
    ir.BPEnergy.set_energy_table(tab)

class RNAConstraintNetwork:
    def __init__(self, seqlen, structures, weights, gcweight):
        self.seqlen = seqlen
        self.structures = list(structures)
        self.weights = list(weights)
        self.gcweight = gcweight
        
    def generate_dependencies(self):
        self.dependencies = list()
                
        bps = set()
        for s in self.structures:
            for bp in s:
                bps.add(bp)
        bps = list(bps)

        self.dependencies = [ [i,j] for (i,j) in bps ]

class RNAConstraintNetworkBasePair(RNAConstraintNetwork):

    def __init__(self, seqlen, structures, weights, gcweight):
        super().__init__(seqlen, structures, weights, gcweight)
        self.generate_cn_basepair_model()

    ## generate constraint network for the base pair model
    def generate_cn_basepair_model(self):
        self.generate_dependencies()

        # generate constraints and functions; assign them to bags
        self.constraints = [ ([i,j], [ir.ComplConstraint(i,j)]) for [i,j] in self.dependencies ]
        self.functions = list()
        
        for k,structure in enumerate(self.structures):
            self.functions.extend( [ ( [i,j], [ir.BPEnergy(i,j,self.weights[k])] ) for (i,j) in structure ] )

        gc_control_funs = [ ( [i], [ir.GCControl(i,self.gcweight)] ) for i in range(self.seqlen) ]
        self.functions.extend( gc_control_funs )
        
class RNATreeDecomposition:
    def __init__(self, cn, *, strategy):
        self.cn = cn

        # from dependencies generate list of binary edges
        bindependencies  = self.expand_to_cliques(cn.dependencies)

        # generate tree decomposition -> bags, edges  (beware: translate between 0/1-based)
        bags,edges = treedecomp.makeTD(cn.seqlen, bindependencies, strategy = strategy)

        self.edges = edges
        self.bags = list(map(set,bags))

    @staticmethod
    def expand_to_cliques(dependencies):
        bindeps = list()
        for d in dependencies:
            bindeps.extend( itertools.combinations(d,2) )
        return bindeps

    ## get (first) index of bag that contains all variables
    def find_bag(self,bvars):
        return next( ( i for i,bag in enumerate(self.bags) if all( x in bag for x in bvars ) ), None )

    # asisgn constraints or functions to bags
    # @returns list where constraints are placed at corresponding bag indices
    def assign_to_bags(self,constraints):
        bagconstraints = { i:[]  for i in range(len(self.bags)) }
        for (cvars,ccons) in constraints:
            bagconstraints[self.find_bag(cvars)].extend(ccons)
        return bagconstraints

    @staticmethod
    def toposort(n,adj):
        visited = set()
        sorted = list()

        def toposort_helper(i):
            visited.add(i)
            for j in adj[i]:
                if not j in visited:
                    toposort_helper(j)
            sorted.append(i)

        for i in range(n):
            if not i in visited:
                toposort_helper(i)
        return sorted[::-1]

    def construct_cluster_tree(self):
        bagconstraints = self.assign_to_bags(self.cn.constraints)
        bagfunctions = self.assign_to_bags(self.cn.functions)

        ct = ir.ClusterTree(self.cn.seqlen, 4);

        ## perform topological sort
        adj = { i:[] for i in range(len(self.bags))}
        for [i,j] in self.edges:
            adj[i].append(j)

        sorted_bags = self.toposort(len(self.bags),adj)

        children = set() # keep record of all non-root nodes (which have been seen as children)
        for bagidx in sorted_bags:
            if not bagidx in children:
                # enumerate subtree
                stack = [(None,bagidx)]
                while stack:
                    (p,i) = stack[-1]
                    stack = stack[:-1]
                    if p==None:
                        cluster = ct.add_root_cluster(list(self.bags[i]))
                    else:
                        cluster = ct.add_child_cluster(p,list(self.bags[i]))

                    for xcon in bagconstraints[i]:
                        ct.add_constraint(cluster, xcon)
                    for xfun in bagfunctions[i]:
                        ct.add_function(cluster, xfun)

                    #print(self.bags[i],list(map(lambda x:x.vars(),bagconstraints[i])))

                    for j in adj[i]:
                        children.add(j)
                        stack.append((cluster,j))
        return ct
