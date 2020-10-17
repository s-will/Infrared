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

from libinfrared import ComplConstraint, BPEnergy, StackEnergy, GCControl, SameComplClassConstraint, DifferentComplClassConstraint
import libinfrared as libir


## @brief exception to signal inconsistency, when consistency would be required
class ParseError(RuntimeError):
    def __init__(self, arg):
        self.args = [arg]


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
                if len(stack[op]) == 0:
                    raise ParseError("Unbalanced RNA dot-bracket structure reading "+cl+".")
                j = stack[op].pop()
                bps[i] = j
                bps[j] = i

    for op in opening:
        if len(stack[op]) > 0:
            raise ParseError("Unbalanced RNA dot-bracket structure reading "+op+".")

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

    if any( len(s) != len(structures[0]) for s in structures[1:] ):
        raise IOError("Read structures of unequal length")

    return structures


# ------------------------------------------------------------
# further RNA specific definitions

## @brief Convert integer (variable value) to nucleotide
## @note encoding A=0, C=1, G=2, U=3
def value_to_nucleotide(x):
    return "ACGU"[x]

## @brief Convert list of integers (variable values) to string
## (sequence) of nucleotides
def values_to_sequence(xs):
    return "".join(map(value_to_nucleotide, xs))

## Parameters for the base pair model (magic params from the Redprint paper)
def_params_bp = { "GC_IN": -2.10208, "AU_IN": -0.52309, "GU_IN": -0.88474,
                  "GC_TERM": -0.09070, "AU_TERM": 1.26630, "GU_TERM": 0.78566 }

## Parameters for the stacking model (magic params from the Redprint paper)
def_params_stacking = { "AUAU": -0.18826, "AUCG": -1.13291, "AUGC": -1.09787,
                        "AUGU": -0.38606, "AUUA": -0.26510, "AUUG": -0.62086,
                        "CGAU": -1.11752, "CGCG": -2.23740, "CGGC": -1.89434,
                        "CGGU": -1.22942, "CGUA": -1.10548, "CGUG": -1.44085,
                        "GUAU": -0.55066, "GUCG": -1.26209, "GUGC": -1.58478,
                        "GUGU": -0.72185, "GUUA": -0.49625, "GUUG": -0.68876 }

## @brief set the bp energy table for Infrared
##
## @param params dictionary of parameters or a table of the parameters
## as expected by infrared::rnadesign
def set_bpenergy_table(params = def_params_bp):
    if type(params) == dict:
        params = list(map(lambda x: params[x],
                          [ "AU_IN", "GC_IN", "GU_IN",
                            "AU_TERM", "GC_TERM", "GU_TERM" ] ))
    libir.BPEnergy.set_energy_table(params)

## @brief set the stacking energy table for Infrared
##
## @param params dictionary of parameters or a table of the parameters
## as expected by infrared::rnadesign
def set_stacking_energy_table(params = def_params_stacking):
    if type(params) == dict:
        params = list(map(lambda x: params[x],
                          [ "AUAU", "AUUA",
                            "AUCG", "AUGC",
                            "AUGU", "AUUG",

                            "CGAU", "CGUA",
                            "CGCG", "CGGC",
                            "CGGU", "CGUG",

                            "GUAU", "GUUA",
                            "GUCG", "GUGC",
                            "GUGU", "GUUG" ] ))
    libir.StackEnergy.set_energy_table(params)


if __name__ == "__main__":
    pass
