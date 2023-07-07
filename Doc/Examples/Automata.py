# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Automata
# ## Forbid/enforce sequence motifs
#
# We showcase an Infrared feature network model that forbids (or enforces) the occurrence of any (respectively one) of several given input motifs in randomly generated sequences.
# The Infrared model is be based on a finite automaton that we construct from the input motifs/words.
#
# For example, let's generate sequences that avoid any of three restriction sites of EcoRI, BamHI, and HindIII

restriction_sites = ["GAATTC", "GGATCC", "AAGCTT"]

import infrared as ir
from collections import defaultdict


# A *Deterministic Finite Automaton (DFA)* is a $5$-tuple $(\mathcal{Q},\Sigma,\delta,q_0,\mathcal{Q}_F)$ with
#
# * Finite set of states $\mathcal{Q}$;
# * Finite set of symbols $\Sigma$;
# * Transition function $\delta:\mathcal{Q}\times\Sigma\to\mathcal{Q}$;
# * Initial state $q_0\in\mathcal{Q}$;
# * Set of final, accepting states $\mathcal{Q}_F\subset\mathcal{Q}$.
#

# ## Construction of an Aho-Corasick automaton
#
# We construct a specific DFA, called *Aho-Corasick automaton*, that accepts strings that contain one of several given words. The general idea of this construction is that the automaton contains one state for every true word prefix; in addition there is a terminal state $. Then, we introduce transitions from states so that the automaton always remembers the longest infix up to current position that can still be extended to an accepted word.

# +
def maxoverlap(s,t):
    for i in reversed(range(min(len(s),(len(t)))+1)):
        if s[-i:] == t[:i]:
            return i
    return 0

def construct_automaton(words):
    #
    Sigma = set()
    for w in words:
        Sigma = Sigma.union(w)

    states = set()
    for w in words:
        for i in range(0,len(w)):
            states.add(w[:i])

    states.add("$")

    # construct transition table as dictionary delta[state][c] = nextstate
    delta = defaultdict(dict)

    for state in states:
        if state=="$":
            for c in Sigma:
                delta[state][c]="$"
            continue
        for c in Sigma:
            s = state+c

            if s in words:
                delta[state][c] = "$"
                continue

            if s in states:
                delta[state][c] = s
                continue

            # find smallest state with largest overlap
            xs = []
            for state2 in states:
                k = maxoverlap(s,state2)
                xs.append((state2,k))
            if len(xs)>0:
                maxk = max(x[1] for x in xs)
                maxstates = [x[0] for x in xs if x[1]==maxk]
                state2=min(maxstates, key=lambda x: len(x))
                delta[state][c]=state2
    return (states,delta)


# -

# ### Construct for specific instance and visualize
#
# As preparation step, to finally forbid the restriction sites, we are first construting a DFA that recognizes them. 

automaton = construct_automaton(restriction_sites)

# +
import graphviz

def draw_automaton(automaton):
    states,delta = automaton
    
    G=graphviz.Digraph(engine="dot",
    graph_attr=dict(),
    node_attr=dict(fontsize="22pt", fontname="Helvetica",
                   penwidth="3.0", fontcolor="black",
                   shape="rectangle", style="rounded",
                   height="0.4", margin="0.03"),
    edge_attr=dict(color="gray30", penwidth="1.5", arrowsize="0.8"))
    
    for state in states:
        G.node(state)
    for src in delta:
        tgts = defaultdict(list)
        for c in delta[src]:
            tgt = delta[src][c]
            tgts[tgt].append(c)
        for tgt,c in tgts.items():
            G.edge(src,tgt,label=','.join(sorted(c)))
    
    return G


# -

draw_automaton(automaton)

# ## Automaton feature network
#
# We 'implement' the automaton as Infrared model, by modeling the state at each position by a variable (whose value is the state index) and the transitions by constraints.
# For use in our model, we define a constraint to express the transition function $\delta$. To express our model in Infrared, states (and characters) have to be encoded as integers. Thus, before we can express the model, we translate our automaton (including its transition function) to this encoding. 

# +
ir.def_constraint_class('Trans',
    lambda i,delta,var: var([('X',i),('Y',i),('Y',i+1)]),
    lambda xi,yi,yip1,delta: yip1 == delta[yi][xi])

def automaton_to_int_encoding(automaton, Sigma):
    """Convert automaton to an integer encoding
    
    Returns:
        (revintstates, intdelta), where
          * revinstates represents the states as a mapping from state labels to integers
          * intdelta is the delta function on integers for states and characters (w.r.t. Sigma)
          
    """
    states, delta = automaton
    
    states=list(states)
    termidx = states.index('$')
    states[termidx],states[-1] = states[-1],states[termidx]
    
    numstates = len(states)
    
    revintstates = {state:x for x,state in enumerate(states)}
        
    intdelta = defaultdict(dict)
    for src in delta:
        for c in delta[src]:
            tgt = delta[src][c]
            intdelta[revintstates[src]][Sigma.index(c)] = revintstates[tgt]
    
    return revintstates, intdelta

    
def automaton_model(n, Sigma, automaton, model=None, negated=False):
    """
    Construct automaton model
    Args:
       Sigma: alphabet
       automaton: actually, only the pair of states and delta of an AC automaton, implicit start at ''
       negated: if False, choose '$' as accepting; if True, every other state
       model: if this is not None, extend the given model by the automaton.
           Given models must define n variables X_i, which encode alphabet characters.
    Returns:
        The constructed model
    """
    
    # sort the alphabet
    Sigma = sorted(Sigma)
    
    # unless a model is already provided, construct it
    # ()
    if model is None:
        model = ir.Model(n,len(Sigma))
    
    states, delta = automaton_to_int_encoding(automaton, Sigma)
    
    q0 = states['']
    qf = states['$']
    
    numstates = len(states)
    
    model.add_variables(1,   (q0,q0),'Y')
    model.add_variables(n-1, numstates ,'Y')
    if not negated:
        model.add_variables(1,   (qf,qf), 'Y')
    else:
        model.add_variables(1,   (0,qf-1), 'Y')
    
    model.add_constraints( Trans(i,delta,model.idx) for i in range(n) )
    
    return model


# -

# ## Model for our example
#
# Construct the model for our example input words and use it to generate (uniformly drawn) sequences that avoid the input words!

# +
Sigma = "ACGT"
n=100
def a_to_seq(assignment, n, Sigma):
    return "".join([Sigma[i] for i in sample.values()[:n]])

model = automaton_model(n, Sigma, automaton, negated=True)
# -

# generate sequences and test whether they contain any input word (if above negated=True, then they must not occur; otherwise one of them must occur)

# +
sampler = ir.Sampler(model)

for _ in range(10):
    sample = sampler.sample()
    seq = a_to_seq(sample,n,Sigma)
    print(seq)
    ## Test whether input words are contained
    for w in restriction_sites:
        if w in seq: 
            print(f"    {w} occurs")
