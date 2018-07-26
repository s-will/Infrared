#!/usr/bin/env python3

import itertools
import infrared as ir

import treedecomp
import rna_support as rna

def val2nucl(x):
    return "ACGU"[x]
def values2seq(xs):
    return "".join(map(val2nucl,xs))

def seed(x):
    ir.seed(x)

def set_bpenergy_table(tab):
    ir.BPEnergy.set_energy_table(tab)


def accumulate_dict(xys):
    d = { x:[] for x,y in xys }
    for (x,y) in xys:
        d[x].append(y)
    return d

class RNAConstraintNetwork:
    def __init__(self, seqlen, structures, weights, gcweight):
        self.seqlen = seqlen
        self.structures = list(structures)
        self.weights = list(weights)
        self.gcweight = gcweight
        
class RNAConstraintNetworkBasePair(RNAConstraintNetwork):

    def __init__(self, seqlen, structures, weights, gcweight):

        super().__init__(seqlen, structures, weights, gcweight)

        self.generate_cn_basepair_model()
        
        self.compute_compl_classes()
        
    ## generate the dependencies due to base pairs
    ## (makes dependencies unique)
    def generate_bpdependencies(self):
        bps = set()
        for s in self.structures:
            for bp in s:
                bps.add(bp)
        bps = list(bps)

        return [ [i,j] for (i,j) in bps ]

    ## generate constraint network for the base pair model
    def generate_cn_basepair_model(self):
        self.dependencies = self.generate_bpdependencies()

        # generate constraints and functions; assign them to bags
        self.constraints = [ ([i,j], [ir.ComplConstraint(i,j)]) for [i,j] in self.dependencies ]
        self.functions = list()
        
        for k,structure in enumerate(self.structures):
            self.functions.extend( [ ( [i,j], [ir.BPEnergy(i,j,self.weights[k])] ) for (i,j) in structure ] )

        gc_control_funs = [ ( [i], [ir.GCControl(i,self.gcweight)] ) for i in range(self.seqlen) ]
        self.functions.extend( gc_control_funs )


    ## compute complementarity classes for each connected component,
    ## i.e. the components of the bipartitiion induced by the
    ## complementarity constraints on base pairs
    ##
    ## @writes result in self.compl_classes
    def compute_compl_classes(self):
        ## perform depth-first traversal
        self.compl_classes = dict()

        visited = set()

        other_ends = accumulate_dict(self.dependencies + list( map(list,map(reversed,self.dependencies)) ))

        color = 1
        
        for x in range(self.seqlen):
            if x in visited: continue
            if x not in other_ends: continue # if there is no base pair

            ## new component, color it
            stack = [x]
            self.compl_classes[x] = color
            color+=1

            while stack:
                x = stack.pop()
                if x in visited: continue
                visited.add(x)
                for y in other_ends[x]:
                    self.compl_classes[y]= -self.compl_classes[x]
                    stack.append(y)

            #print(self.compl_classes)

    # def gen_redundant_bpconstraints(self,vars):
    #     if not hasattr(self,"compl_classes"): return []
    #     constraints = list()
    #     for j in range(1,len(vars)):
    #         # find component
    #         for cc in self.compl_classes:
    #             if vars[j] in cc:
    #                 for i in range(0,j):
    #                     if vars[j] in cc:
    #                         if cc[vars[i]]==cc[vars[j]]:
    #                             #same component
    #                             constraints.append(ir.SameComplClassConstraint(i,j))
    #                             break
    #                         else:
    #                             #different classes
    #                             constraints.append(ir.DifferentComplClassConstraint(i,j))
    #                             break
    #     return constraints

class RNATreeDecomposition:
    def __init__(self, cn, *, add_redundant_constraints=True, strategy):
        self.cn = cn

        self.add_redundant_constraints = add_redundant_constraints

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
    def find_all_bags(self,bvars):
        return [ i for i,bag in enumerate(self.bags) if all( x in bag for x in bvars ) ]

    def find_bag(self,bvars):
        bags = self.find_all_bags(bvars)
        if bags is not []:
            return bags[0]
        else:
            return None

    # asisgn constraints or functions to bags
    # @returns list where constraints are placed at corresponding bag indices
    def assign_to_bags(self,constraints):
        bagconstraints = { i:[]  for i in range(len(self.bags)) }
        for (cvars,ccons) in constraints:
            bagconstraints[self.find_bag(cvars)].extend(ccons)
        return bagconstraints

    def assign_to_all_bags(self,constraints):
        bagconstraints = { i:[]  for i in range(len(self.bags)) }
        for (cvars,ccons) in constraints:
            for bidx in self.find_all_bags(cvars):
                bagconstraints[bidx].extend(ccons)
        return bagconstraints

    ## Fills in constraints due to complementarity to all bags. Adds
    ## SameComplClass or DifferentComplClass constraints, wherever a
    ## variable would have to be enumerated unnconstrained otherwise
    ## @param classes complementarity classes
    ## @param existing_constraints
    ## @param[in,out] bagconstraints
    ##
    ## avoids to add constraints which already exist
    def fillin_class_constraints(self, classes, existing_constraints, bagconstraints):
        existing = set( (vars[0],vars[1]) for (vars,constr) in existing_constraints )
        for bagidx,bag in enumerate(self.bags):
            bagvars=sorted(list(bag))
            for j in bagvars[1:]:
                if all( (i,j) not in existing for i in range(0,j) ):
                    i = bagvars[0]
                    if classes[i] == classes[j]:
                        #print("Add same:",i,j)
                        bagconstraints[bagidx].append(ir.SameComplClassConstraint(i,j))
                    elif classes[i] == - classes[j]:
                        #print("Add diff:",i,j)
                        bagconstraints[bagidx].append(ir.DifferentComplClassConstraint(i,j))

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
        if self.add_redundant_constraints:
            bagconstraints = self.assign_to_all_bags(self.cn.constraints)
            self.fillin_class_constraints(self.cn.compl_classes,
                                          self.cn.constraints,
                                          bagconstraints)
        else:
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
                    bagvars = sorted(list(self.bags[i]))

                    if p==None:
                        cluster = ct.add_root_cluster(bagvars)
                    else:
                        cluster = ct.add_child_cluster(p,bagvars)

                    for xcon in bagconstraints[i]: #  + self.cn.gen_redundant_bpconstraints(bagvars)
                        ct.add_constraint(cluster, xcon)

                    for xfun in bagfunctions[i]:
                        ct.add_function(cluster, xfun)

                    #print(self.bags[i],list(map(lambda x:x.vars(),bagconstraints[i])))

                    for j in adj[i]:
                        children.add(j)
                        stack.append((cluster,j))
        return ct
