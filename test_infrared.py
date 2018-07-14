#!/usr/bin/env python3

import infrared as ir
from collections import Counter

def val2nucl(x):
    return "ACGU"[x]
def values2seq(xs):
    return "".join(map(val2nucl,xs))

class SoftComplConstr(ir.Function):
    def __init__(self,i,j,w):
        ir.Function.__init__(self,[i,j])
        self.w = w
        self.i = i
        self.j = j
    def __call__(self,a):
        v = a.values()
        if v[self.i] == 3-v[self.j]:
            e = -1
        else:
            e = 0
        return pow(self.w, -e)

w_C = 5
w_CG = 2
w_E = 10

# a = ir.Assignment(10)

seqlen=30

ct = ir.ClusterTree(seqlen,4);

bpairs = [(1,9), (2,8), (3,5), (6,7), (11,13), (0,14)]
structure= [ '.' for i in range(0,seqlen) ]
for (i,j) in bpairs:
    structure[i]='('
    structure[j]=')'
structure = "".join(structure)

for (i,j) in bpairs:
    cluster = ct.add_root_cluster([i,j])
    #ct.add_constraint( cluster, ir.ComplConstraint(i,j) ) 
    ct.add_function( cluster, SoftComplConstr(i,j, w_C) ) 
    #ct.add_function( cluster, ir.BPEnergy(i,j, w_E) ) 
    ct.add_function( cluster, ir.CGControl(i, w_CG) )
    ct.add_function( cluster, ir.CGControl(j, w_CG) )

paired = set( [ i for (i,j) in bpairs ] + [ j for (i,j) in bpairs ] )

for i in range(0,seqlen):
    if i in paired: continue
    cluster = ct.add_root_cluster([i])
    ct.add_function( cluster, ir.CGControl(i, w_CG) )

ct.evaluate()

## statistics
counters=[Counter() for i in range(0,seqlen)]

print(structure)
for x in range(0,20):
    sample = ct.sample()
    seq = values2seq(sample.values())
    print(seq)

    ## statistics
    for i in range(0,seqlen):
        counters[i][seq[i]] += 1


print(sum(counters[1:],counters[0]))



# class SimpleBPEnergy(ir.Function):
#     def __init__(self,i,j):
#         self.i = i
#         self.j = j

#     def eval(self, a):
#         if a[i] == 3-a[j]: return 0
#         else: return 1e-4

