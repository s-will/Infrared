# This file is part of InfraRed source code.

## @file
##
## @brief Code related to automaton integration in InfraRed

from automata.fa.nfa import NFA
from automata.fa.dfa import DFA

from .infrared import def_constraint_class

SYMBOLS = {'0', '1', '2', '3'}
NUC_TO_SYM = {'A':'0', 'a':'0', 'C':'1', 'c':'1', 'G':'2', 'g':'2', 'U':'3', 'u':'3'}

#TODO: complete
## @brief Return nfa given a regular expression
def regex_to_nfa(regex):
    pass

def def_state(i,j):
    return '{}-{}'.format(i,j)

## @brief DFA recognizes one of input words
def words_to_dfa(words):
    states = {'s', 'f'}
    transitions = {'s':{c: {'s'} for c in SYMBOLS}, 'f':{c:'f' for c in SYMBOLS}}
    for ind, word in enumerate(words):
        for i, c in enumerate(word):
            c = NUC_TO_SYM[c]
            states.add(def_state(ind, i))
            if i == 0:
                transitions['s'][c].add(def_state(ind,i))
            elif i == len(word)-1:
                transitions.setdefault(def_state(ind,i-1),{}).setdefault(c, set()).add('f')
            else:
                transitions.setdefault(def_state(ind,i-1),{}).setdefault(c, set()).add(def_state(ind,i))
    nfa = NFA(states=states, input_symbols=SYMBOLS, transitions=transitions, initial_state='s', final_states={'f'})
    dfa = DFA.from_nfa(nfa)
    return rename_dfa(dfa.minify())

## @Rename the state of given dfa
def rename_dfa(dfa):
    states_map = {s:'q'+str(i) for i,s in enumerate(dfa.states)}
    states = {states_map[s] for s in dfa.states}
    initial_state = states_map[dfa.initial_state]
    finial_states = {states_map[s] for s in dfa.final_states}
    transitions = {states_map[s]: {c:states_map[q] for c, q in v.items()} for s,v in dfa.transitions.items()}
    return DFA(states=states, input_symbols=SYMBOLS, transitions=transitions, initial_state=initial_state, final_states=finial_states)

## @brief Return constraints given a dfa
##
## States Id are name
def dfa_to_constraints(dfa, start, end, state_start):
    states_nb = len(dfa.states)

    def to_dfa_state(state):
        return 'q'+str(state)

    # Transition
    def_constraint_class('Transition', lambda i: [start+i, state_start+i, state_start+i+1], lambda x, s1, s2: rule_stop_at_final(str(x), to_dfa_state(s1), to_dfa_state(s2), dfa), module=__name__)

    # Start State
    def_constraint_class('StartState', lambda i: [i], lambda s: to_dfa_state(s) == dfa.initial_state, module= __name__)
    # Final State
    def_constraint_class('FinalState', lambda i: [i], lambda s: to_dfa_state(s) in dfa.final_states, module=__name__)

    return [states_nb]*(end-start+2), [Transition(i= i-start) for i in range(start, end+1)] + [StartState(i=state_start), FinalState(i=state_start+end-start)]


## @brief Stay at final states once arrived
##
##
def rule_stop_at_final(x, s1, s2, dfa):
    if s1 in dfa.final_states:
        return s1 == s2
    else:
        return dfa.transitions[s1][x] == s2


def words_to_accept(words, end, state_start, start=0):
    dfa = words_to_dfa(words)
    variables, constraints = dfa_to_constraints(dfa, start, end, state_start)
    return variables, constraints

if __name__ == '__main__':
    dfa = DFA(
    states={'q0', 'q1', 'q2'},
    input_symbols={'0', '1'},
    transitions={
        'q0': {'0': 'q0', '1': 'q1'},
        'q1': {'0': 'q0', '1': 'q2'},
        'q2': {'0': 'q2', '1': 'q1'}
    },
    initial_state='q0',
    final_states={'q1'}
    )
