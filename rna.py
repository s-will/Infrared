#!/usr/bin/env python3

# some RNA related functions

import re

def read_inp(inpfh):
    """read structures from inp file
    @param fh input file handle to inp file
    """

    structures=list()

    for line in inpfh:
        if re.match(re.compile("[\(\)\.]"), line, flags=0):
            structures.append(line.rstrip('\n'))
        elif re.match(";", line, flags=0):
            break

    return structures

def parseRNAStructure(structure):
    """Parse RNA structure including pseudoknots
    @param structure
    @returns array
    """

    opening = "([{<"
    closing = ")]}>"

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

def parseRNAStructureBps(structure):
    """Parse dot bracket structure and return list of base pairs"""
    s = parseRNAStructure(structure)
    bps = list()
    for i,j in enumerate(s):
        if i != -1 and i<j:
            bps.append( (i,j) )
    return bps


def unique_edges(xs):
    d = { (min(x,y),max(x,y)):1 for (x,y) in xs }
    return list(d.keys())

def structure_to_basepair_dependencies(structure,edges=[]):
    """Convert structure to list of edges in basepair model; node indices are 0-based!"""

    for (i,j) in enumerate(structure):
        if i<j:
            edges.append((i,j))
    

def structure_to_stacking_dependencies(structure,edges=[]):
    """Convert structure to list of edges in stacking model; node indices are 0-based!"""

    for (i,j) in enumerate(structure):
        if i+1<j-1 and structure[i+1]==j-1:
            edges.extend([(i,j),(i+1,j-1),(i+1,j),(i,j-1),(i,i+1),(j,j-1)])
    
if __name__ == "__main__":
    pass
