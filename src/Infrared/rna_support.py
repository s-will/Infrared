#!/usr/bin/env python3

#
# InfraRed ---  A generic engine for Boltzmann sampling over constraint networks
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

## @file
##
## @brief Some RNA related functions
##
## Loose library for some common tasks for RNA specific code
##

import re
import collections

from libinfrared import ComplConstraint,BPEnergy,GCControl


## @brief Parse RNA structure including pseudoknots
##
## @param structure dot bracket string of RNA structure
## @param opening specifies the recognized opening bracket symbols
## @param closing specifies closing bracket symbols (corresponding to opening)
##
## @note Characters not in opening or closing are considered unpaired.
##
## @return array bps encoding the structure like bps[i]=j for each bp
## {i,j}
def parseRNAStructure(structure, *, opening = "([{<", closing = ")]}>"):
    stack = { op:list() for op in opening }
    bps = [-1]*len(structure)

    for i,c in enumerate(structure):
        for (op,cl) in zip(opening,closing):
            if c==op:
                stack[op].append(i)
            elif c==cl:
                j = stack[op].pop()
                bps[i] = j
                bps[j] = i

    return bps

## @brief Parse RNA structure, returning list of base pairs
##
## @param structure dot bracket string of RNA structure
## @see parseRNAStructure
##
## @returns list of base pairs (i,j)
def parseRNAStructureBps(structure, **kwargs):
    s = parseRNAStructure(structure, **kwargs)
    bps = list()
    for i,j in enumerate(s):
        if i != -1 and i<j:
            bps.append( (i,j) )
    return bps

## @brief Check complementarity of nucleotides
## @param x nucelotide
## @param y nucelotide
## @return whether complementary (A-U,C-G, or G-U)
def is_complementary(x,y):
    compls=set(["AU","CG","GU"])
    return x+y in compls or y+x in compls

## @brief Get invalid base pairs
##
## @param seq sequence string
## @param struc structure as dot bracket string or base pair list
##
## @return list of base pairs that violate the complementarity constraints
def invalid_bps(seq, struc):
    if type(struc)==str:
        bps = parseRNAStructureBps(struc)
    else:
        bps = struc
    invalids = []
    for (i,j) in bps:
        if not is_complementary(seq[i],seq[j]):
            invalids.append((i,j))

    return invalids

## @brief Check whether sequence satisfies complementarity constraints due to structure
## @param seq sequence string
## @param struc structure as dot bracket string or base pair list
## @return whether valid
def is_valid(seq, struc):
    return invalid_bps(seq,struc) == []


## @brief Convert structure to list of edges in basepair model
##
## @param structure array representation of RNA structure
## @param[in,out] edges list of edges
##
## Dependencies are appendes to edges
##
## @note node indices are 0-based
def structure_to_basepair_dependencies(structure,edges=[]):
    for (i,j) in enumerate(structure):
        if i<j:
            edges.append((i,j))

## @brief Convert structure to list of edges in stacking model
##
## @param structure array representation of RNA structure
## @param[in,out] edges list of edges
##
## Dependencies are appendes to edges
## @see structure_to_basepair_dependencies
def structure_to_stacking_dependencies(structure,edges=[]):
    for (i,j) in enumerate(structure):
        if i+1<j-1 and structure[i+1]==j-1:
            edges.extend([(i,j),(i+1,j-1),(i+1,j),(i,j-1),(i,i+1),(j,j-1)])

## @brief Compute GC content in a sequence
## @param seq sequence string
## @return gc content (as ratio in [0,1])
def GC_content(seq):
    c = collections.Counter(seq)

    gc = 0
    for x in ['C','G']:
        if x in c:
            gc += c[x]

    gc = gc  / len(seq)

    return gc

## @brief Make edges unique (considering symmetry)
## @param xs list of edges
## @return unique list of edges
##
## keeps only one of (x,y) and (y,x), namely (x,y) if x<y
def unique_edges(xs):
    d = { (min(x,y),max(x,y)):1 for (x,y) in xs }
    return list(d.keys())


## @brief Read multiple structures from file in inp format
##
## @param inpfh input file handle
##
## @return list of the structures as dot bracket strings
##
## inp format is a file format to specify instances for multi-target
## RNA design, e.g. used in the benchmarks of Modena and RNAblueprint
def read_inp(inpfh):
    structures=list()

    for line in inpfh:
        if re.match(re.compile("[\(\)\.]"), line, flags=0):
            structures.append(line.rstrip('\n'))
        elif re.match(";", line, flags=0):
            break

    return structures



if __name__ == "__main__":
    pass
