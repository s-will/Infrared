#!/usr/bin/env python3

#
# InfraRed -  A generic engine for Boltzmann sampling over constraint networks
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

## @file rna.py
#  @brief Functionality for RNA-related Infrared applications
#
#  @code
#    import infrared.rna as rna
#  @endcode
#

##
#  @package infrared.rna
#  @copydoc rna.py


import re
import collections
import math

from . import libinfrared
from .infrared import def_constraint_class, def_function_class, WeightedFunction

# @brief exception to signal inconsistency, when consistency would be required
class ParseError(RuntimeError):
    def __init__(self, arg):
        self.args = [arg]

# #####
# define some constraints and functions for RNA design
#
# Note that these classes are solely defined by calling the functions
# def_constraint_class and def_function_class.
#
# Here, we define the classes with comments only for the purpose of getting
# clean doxygen documentation. These classes will be overwritten by the
# def_*_class functions.
#

class BPComp(libinfrared.Constraint):
    """
    Constrain complementarity of the base pair (i,j)

    ```
    BPComp(i,j)
    ```
    The constraint is satisfied if values at positions (i,j) form a valid canonical base pair, _i.e._ {(A,U), (C,G), (G,U)}.
    """
## @cond PRIVATE
_bpcomp_tab = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2)]
def_constraint_class('BPComp', lambda i, j: [i, j],
                     lambda x, y: (x, y) in _bpcomp_tab,
                     module=__name__)
## @endcond



class NotBPComp(libinfrared.Constraint):
    """Constraint for negation of BPComp

    ```
    NotBPComp(i,j)
    ```
    The constraint is satisfied if values at positions (i,j) DO NOT form a valid canonical base pair.

    @see BPComp
    """
## @cond PRIVATE
def_constraint_class('NotBPComp', lambda i, j: [i, j],
                     lambda x, y: (x, y) not in _bpcomp_tab,
                     module=__name__)
## @endcond

class GCCont(WeightedFunction):
    """
    Function for (position-wise)  GC content

    ```
    GCCont(i)
    ```
    GCCont is an Infrared Function to count GCCont at position i, 1 if the value is C or G, 0 otherwise.
    """

## @cond PRIVATE
def_function_class('GCCont', lambda i: [i],
                   lambda x: 1 if x == 1 or x == 2 else 0,
                   module=__name__)
## @endcond

class BPEnergy(WeightedFunction):
    """
    Function for (basepair-wise) BasePair Energy model

    BPEnergy is an Infrared Function to capture BasePair Energy model at base pair (i,j).
    It further takes a boolean value indicating whether the base pair (i,j) is terminal, _i.e._ (i-1,j+1) is not a base pair

    ```
    bps = parse(target)
    bpFunctions = [BPEnergy(i,j, (i-1,j+1) not in bps) for (i,j) in bps]
    ```
    """

## @cond PRIVATE
def_function_class('BPEnergy', lambda i, j, is_terminal: [i, j],
                   lambda x, y, is_terminal: _bpenergy(x, y, is_terminal),
                   module=__name__)
## @endcond

class StackEnergy(WeightedFunction):
    """
    Function for Stack Energy model

    StackEnergy is an Infrared Function to capture Stack Energy model for a base pair stack (i, j, i+1, j-1).

    ```
    bps = parse(target)
    stackFunctions = [StackEnergy(i,j) for (i,j) in bps if (i+1,j-1) in bps]
    ```
    """

## @cond PRIVATE
def_function_class('StackEnergy', lambda i, j: [i, j, i+1, j-1],
                   lambda x, y, x1, y1: _stackenergy(x, y, x1, y1),
                   module=__name__)
## @endcond

class SameComplClassConstraint(libinfrared.Constraint):
    """
    Constrain nucleotides to be in the same complementarity class

    ```
    SameComplClassConstraint(i,j)
    ```
    The constraint is satisfied if values at positions (i,j) are from the same complenentarity class, i.e. either both in {A,G} or both in {C,U}.
    """

## @cond PRIVATE
def_constraint_class('SameComplClassConstraint', lambda i, j: [i, j],
                     lambda x, y: x & 1 == y & 1,
                     module=__name__)
## @endcond

class DifferentComplClassConstraint(libinfrared.Constraint):
    """
    Constrain nucleotides to be in different complementarity classes

    ```
    DifferentComplClassConstraint(i,j)
    ```
    Negation of SameComplClassConstraint.
    The constraint is satisfied if values at positions (i,j) are NOT from the same complenentarity class.
    """

## @cond PRIVATE
def_constraint_class('DifferentComplClassConstraint', lambda i, j: [i, j],
                     lambda x, y: x & 1 != y & 1,
                     module=__name__)
## @endcond

def parse_array(structure, *, opening="([{<", closing=")]}>"):
    """Parse RNA structure including pseudoknots

    Args:
        structure: dot bracket string of RNA structure
        opening: specifies the recognized opening bracket symbols
        closing: specifies closing bracket symbols (corresponding to opening)
    Returns:
        array bps encoding the structure like bps[i]=j for each bp {i,j}

    Notes:
        Characters not in opening or closing are considered unpaired.
    """
    stack = {op: list() for op in opening}
    bps = [-1]*len(structure)

    for i, c in enumerate(structure):
        for (op, cl) in zip(opening, closing):
            if c == op:
                stack[op].append(i)
            elif c == cl:
                if len(stack[op]) == 0:
                    raise ParseError(
                        "Unbalanced RNA dot-bracket structure reading "+cl+".")
                j = stack[op].pop()
                bps[i] = j
                bps[j] = i

    for op in opening:
        if len(stack[op]) > 0:
            raise ParseError(
                "Unbalanced RNA dot-bracket structure reading "+op+".")

    return bps


def parse(structure, **kwargs):
    """Parse RNA structure, returning list of base pairs

    Args:
        structure: dot bracket string of RNA structure

    Returns:
        list of base pairs (i,j)

    See Also:
        parse_array
    """
    s = parse_array(structure, **kwargs)
    bps = list()
    for i, j in enumerate(s):
        if i != -1 and i < j:
            bps.append((i, j))
    return bps


def is_complementary(x, y):
    """Check complementarity of nucleotides
    Args:
        x: nucelotide
        y: nucelotide

    Returns:
        whether complementary (A-U,C-G, or G-U)
    """
    compls = set(["AU", "CG", "GU"])
    return x+y in compls or y+x in compls


def invalid_bps(seq, struc):
    """Get invalid base pairs
    Args:
        seq: sequence string
        struc: structure as dot bracket string or base pair list

    Returns:
        list of base pairs that violate the complementarity constraints
    """
    if type(struc) == str:
        bps = parse(struc)
    else:
        bps = struc
    invalids = []
    for (i, j) in bps:
        if not is_complementary(seq[i], seq[j]):
            invalids.append((i, j))

    return invalids


def is_valid(seq, struc):
    """Check whether sequence satisfies complementarity constraints due
    to structure

    Args:
        seq: sequence string
        struc: structure as dot bracket string or base pair list
    Returns:
        whether valid
    """
    return invalid_bps(seq, struc) == []


def structure_to_basepair_dependencies(structure, edges=[]):
    """Convert structure to list of edges in basepair model.
    Dependencies are appendes to edges

    Args:
        structure: array representation of RNA structure
        edges:  edges list of edges [in,out]

    Note:
        node indices are 0-based
    """
    for (i, j) in enumerate(structure):
        if i < j:
            edges.append((i, j))


def structure_to_stacking_dependencies(structure, edges=[]):
    """Convert structure to list of edges in stacking model.
    Dependencies are appendes to edges

    Args:
        structure: array representation of RNA structure
        edges: list of edges [in,out]

    See Also:
        structure_to_basepair_dependencies
    """
    for (i, j) in enumerate(structure):
        if i+1 < j-1 and structure[i+1] == j-1:
            edges.extend([(i, j), (i+1, j-1), (i+1, j),
                          (i, j-1), (i, i+1), (j, j-1)])


def GC_content(seq):
    """Compute GC content in a sequence
    Args:
        seq: sequence string
    Returns:
        gc content (as ratio in [0,1])
    """
    c = collections.Counter(seq)

    gc = 0
    for x in ['C', 'G']:
        if x in c:
            gc += c[x]

    gc = gc / len(seq)

    return gc


def unique_edges(xs):
    """Make edges unique (considering symmetry).
    keeps only one of (x,y) and (y,x), namely (x,y) if x<y

    Args:
        xs: list of edges
    Returns:
        unique list of edges
    """
    d = {(min(x, y), max(x, y)): 1 for (x, y) in xs}
    return list(d.keys())


def read_inp(inpfh):
    """Read multiple structures from file in inp format.
    inp format is a file format to specify instances for multi-target
    RNA design, e.g. used in the benchmarks of Modena and RNAblueprint

    Args:
        inpfh: input file handle

    Returns:
        list of the structures as dot bracket strings
    """
    structures = list()

    for line in inpfh:
        if re.match(re.compile(r"[\(\)\.]"), line, flags=0):
            structures.append(line.rstrip('\n'))
        elif re.match(";", line, flags=0):
            break

    if any(len(s) != len(structures[0]) for s in structures[1:]):
        raise IOError("Read structures of unequal length")

    return structures


# ------------------------------------------------------------
# further RNA specific definitions

def nucleotide_to_value(x):
    """Convert integer (variable value) to nucleotide

    Args:
        x: nucleotide or gap, A, C, G, U, T, -, or .
    Note:
        encoding A=0, C=1, G=2, U/T=3, -/.=4
    """
    return {'A':0,'C':1,'G':2,'U':3,'T':3,'.':4,'-':4}[x]

def value_to_nucleotide(x):
    """Convert integer (variable value) to nucleotide
    Args:
        x: integer
    Note:
        encoding A=0, C=1, G=2, U=3, -=4
    """
    return "ACGU-"[x]

def values_to_seq(xs):
    """Convert list of integers (variable values) to string (sequence) of nucleotides
    Args:
        xs: list of integers
    """
    return "".join(map(value_to_nucleotide, xs))


def ass_to_seq(ass):
    """Convert assignment to sequence string
    Args:
        ass: assignment
    """
    return values_to_seq(ass.values())

# ------------------------------------------------
#
# support for IUPAC

## @cond PRIVATE
_iupac_nucleotides = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'U',
    'U': 'U',
    'R': 'AG',
    'Y': 'CU',
    'S': 'CG',
    'W': 'AU',
    'K': 'GU',
    'M': 'AC',
    'B': 'CGU',
    'D': 'AGU',
    'H': 'ACU',
    'V': 'ACG',
    'N': 'ACGU',
    '.': '-',
    '-': '-'
}
## @endcond

def iupacvalues(symbol):
    return [ nucleotide_to_value(x) for x in _iupac_nucleotides[symbol] ]

# ------------------------------------------------
# RNA energy models / parameters

## @cond PRIVATE

## Base pair energy parameters
#
# used by _bpenergy and Function BPEnergy
# @see set_bpenergy_table
_params_bp = None

## Stacking energy parameters
#
# used by _stackenergy and Function StackEnergy
# @see set_stacking_energy_table
_params_stacking=None

## @endcond


## Default parameters for the base pair model (magic params from the Redprint paper)
def_params_bp = {"GC_IN": -2.10208, "AU_IN": -0.52309, "GU_IN": -0.88474,
    "GC_TERM": -0.09070, "AU_TERM": 1.26630, "GU_TERM": 0.78566}

## Default parameters for the stacking model (magic params from the Redprint paper)
def_params_stacking = {"AUAU": -0.18826, "AUCG": -1.13291, "AUGC": -1.09787,
    "AUGU": -0.38606, "AUUA": -0.26510, "AUUG": -0.62086,
    "CGAU": -1.11752, "CGCG": -2.23740, "CGGC": -1.89434,
    "CGGU": -1.22942, "CGUA": -1.10548, "CGUG": -1.44085,
    "GUAU": -0.55066, "GUCG": -1.26209, "GUGC": -1.58478,
    "GUGU": -0.72185, "GUUA": -0.49625, "GUUG": -0.68876}


def set_bpenergy_table(params=def_params_bp):
    """Set the bp energy table for Function BPEnergy

    Args:
        params: dictionary of parameters or a table of the parameters
        as expected by _bpenergy
    """
    if type(params) == dict:
        params = list(map(lambda x: params[x],
                          ["AU_IN", "GC_IN", "GU_IN",
                           "AU_TERM", "GC_TERM", "GU_TERM"]))
    global _params_bp
    _params_bp = params


def set_stacking_energy_table(params=def_params_stacking):
    """Set the stacking energy table for Function StackEnergy

    Args:
        params: dictionary of parameters or a table of the parameters
        as expected by _stackenergy
    """
    if type(params) == dict:
        params = list(map(lambda x: params[x],
                          ["AUAU", "AUUA",
                           "AUCG", "AUGC",
                           "AUGU", "AUUG",

                           "CGAU", "CGUA",
                           "CGCG", "CGGC",
                           "CGGU", "CGUG",

                           "GUAU", "GUUA",
                           "GUCG", "GUGC",
                           "GUGU", "GUUG"]))
    global _params_stacking
    _params_stacking = params


# run table initialization
set_bpenergy_table()
set_stacking_energy_table()

## @cond PRIVATE

## lookup table for base pair indices
#
# used by _bpenergy and _stackenergy
_bpindex_tab = [[-1, -1, -1, 0],
                [-1, -1, 2, -1],
                [-1, 3, -1, 4],
                [1, -1, 5, -1]]

def _bpenergy(x, y, is_terminal=False):
    """
    @brief Energy of base pair

    Args:
     x: base in internal 0..3 representation
     y: base in internal 0..3 representation
     is_terminal: flag, True if the base pair is terminating a stem

    Returns:
        energy of the base pair according to current bp energy table

    @see set_bpenergy_table()
    """
    bpidx = _bpindex_tab[x][y]
    return (_params_bp[bpidx//2 + (3 if is_terminal else 0)]
            if bpidx >= 0 else -math.inf)


def _stackenergy(x, y, x1, y1):
    """
    @brief Energy of stack of base pairs

    Args:
     x: base in internal 0..3 representation
     y: base in internal 0..3 representation
     x1: base in internal 0..3 representation, inner stacked base pair
     y1: base in internal 0..3 representation, inner stacked base pair

    Returns:
        energy of the base pair stack according to current stacking energy table

    @see set_stacking_energy_table()
    """
    bpidx = _bpindex_tab[x][y]
    bpidx1 = _bpindex_tab[x1][y1]

    if bpidx < 0 or bpidx1 < 0:
        return -math.inf

    if bpidx & 1:
        bpidx, bpidx1 = bpidx1, bpidx
        bpidx1 ^= 1

    return _params_stacking[6 * (bpidx//2) + bpidx1]

## @endcond

if __name__ == "__main__":
    pass
