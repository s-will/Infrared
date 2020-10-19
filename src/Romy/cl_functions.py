#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 18:34:51 2020

@author: amine-abdeljaoued
"""
import RNA
import numpy as np

#Get a simple pair with the position of base pairs from the list format returned by RNA.plist()
def get_pairs(initial_format): 
    a = str(initial_format)
    first=int(a.split(',')[0].split(':')[1].strip())
    second=int(a.split(',')[1].split(':')[1].strip())

    return (first,second)


def cluster_probability(cluster,structures,n,structs_count,s):
    P_Sk = 1/s #Probability of a sequence, assuming uniform distribution here
    P = 0
    for i in cluster:
        """ index_seq =  i//n #Retrieve the index of initial sequence
        P += ( (structs_count[index_seq][structures[i]])/(sum(structs_count[index_seq].values())) )* P_Sk """

        P += ( (structs_count[structures[i]])/(sum(structs_count.values())) )* P_Sk
    
    return P

#Construct the probability split
def base_pair_proba(cluster, structures, base_pairs,n,structs_count,s):
    probas = {}
    p_cluster = cluster_probability(cluster,structures,n,structs_count, s)
    P_Sk = 1/s #Probability of a sequence, assuming uniform distribution here
    for i in cluster:
        for pair in base_pairs[i]:
            if pair not in probas:
                probas[pair]= (P_Sk* ( (structs_count[structures[i]])/(sum(structs_count.values())) ) / p_cluster)
            else:
                probas[pair]+= (P_Sk* ( (structs_count[structures[i]])/(sum(structs_count.values())) ) / p_cluster)
    
    for i in probas.keys():
        if probas[i]>1.0: 
            probas[i]=1.0 #Avoid having probabilities greater than 1 because of rounding
    return probas

#Other version considering that structures are already present with their given probabilities
def base_pair_proba_v2(cluster, structures, base_pairs,n,s):
    probas={}
    p_cluster=len(cluster)/len(structures)
    P_Sk = 1/s #Probability of a sequence, assuming uniform distribution here
    for i in cluster:
        for pair in base_pairs[i]:
            if pair not in probas:
                probas[pair]=0
                #probas[pair]= (P_Sk* fcs[index_fc].pr_structure(structures[i]) / p_cluster)
            else:
                probas[pair]+=1
                #probas[pair]+= (P_Sk* fcs[index_fc].pr_structure(structures[i]) / p_cluster)
    
    for i in probas.keys():
        probas[i] /= (len(cluster)*p_cluster) 
        if probas[i]>1.0: 
            probas[i]=1.0 #Avoid having probabilities greater than 1 because of rounding
    return probas

def cluster_diversity(cluster_index,probas):
    return sum([p*(1-p) for p in probas[cluster_index].values()])
    
def cluster_ensemble_energy(cluster,structures,T,fc,n,s):
    R = RNA.GASCONST/1000
    P_Sk = 1/s #Probability of a sequence, assuming uniform distribution here
    E = 0
    for i in cluster:
        index_fc =  i//n #Retrieve the index of fold coumpound (the initial sequence)
        E+= np.exp(-fc.eval_structure(structures[i])/ (R*T))* P_Sk
    
    return -R*T* np.log(E)
