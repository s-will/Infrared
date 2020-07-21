#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 17:49:34 2020

@author: amine-abdeljaoued
"""


import RNA
import argparse
import numpy as np
import pandas as pd
#from tabulate import tabulate
from sklearn.cluster import KMeans, SpectralClustering, MiniBatchKMeans
import rna_support as rna

from functions import *



def clustering(sequences,k,n=15):

    s = len(sequences) #Number of sequences in alignment
    N=  n*s #Total number of structures

    fcs = [] #Fold compounds objects
    structs = [] #Structures generated
    base_pairs_list = [] #List of parsed structures: base pair pairs 
    diss_matrix=np.zeros((N,N)) #Dissimilarity matrix of base pair distances
    similarity_matrix = np.zeros((N,N))
    clusters = [[] for _ in range(k)] #Clusters

    average_gc=0
    for seq in sequences:
        #Compute GC content
        average_gc+= rna.GC_content(seq) * 100
        
        # create model details
        md = RNA.md()
        # activate unique multibranch loop decomposition
        md.uniq_ML = 1
        # create fold compound object
        fc = RNA.fold_compound(seq, md)
        # compute MFE[value for value in ls[value for value in lst1 if value in lst2] [value for value in lst1 if value in lst2] [value for value in lst1 if value in lst2] [value for value in lst1 if value in lst2] [value for value in lst1 if value in lst2] [value for value in lst1 if value in lst2] [value for value in lst1 if value in lst2] t1 if value in lst2] 
        (ss, mfe) = fc.mfe()
        # rescale Boltzmann factors according to MFE
        fc.exp_params_rescale(mfe)
        
        # compute partition function to fill DP matrices
        (pp,pf) = fc.pf()
        
        elements=0
        #Backtrace multiple structures non-redundantly
        while(elements< n):
            for struc in fc.pbacktrack(N,RNA.PBACKTRACK_NON_REDUNDANT):
                #print(s,"\n")
                if struc not in structs: 
                    structs.append(struc)
                    elements += 1
                    if elements== n: break
        
        fcs.append(fc)
        
    average_gc /= s

    #Parse Vienna notation into a set of pairs using RNA.plist()
    #And get_pairs defined in functions.py  
    for structure in structs:
        struct_set=set()
        for x in RNA.plist(structure,1):
            struct_set.add(get_pairs(x))
        base_pairs_list.append(struct_set)
        
    #Compute the dissimilarity matrix
    for i in range(N):
        for j in range(N):
            A = base_pairs_list[i]
            B = base_pairs_list[j]
            diss_matrix[i][j] = len(A.symmetric_difference(B))
            similarity_matrix[i][j]=len(A.intersection(B))      

    
    #Clustering with k means
    kmeans = MiniBatchKMeans(n_clusters= k, init='k-means++')
    kmeans.fit(diss_matrix)
    for i in range(N):
        clusters[kmeans.labels_[i]].append(i) #Put the structures in their clusters """ """

    return [clusters,structs, base_pairs_list,fcs,s, kmeans,average_gc]
    
    
"""     #Spectral clustering
    spec_clus = SpectralClustering(n_clusters=k,affinity='precomputed').fit(similarity_matrix) 
    for i in range(N):
        clusters[spec_clus.labels_[i]].append(i) """

    
    


def analyze_clusters(clusters, structs, base_pairs_list, fcs, s, n, sequences,k ,T=310.15, gamma=1):
    
    #MEA structure for each cluster
    reps = [] #Representatives for each cluster
    probas=[]
    for cluster in clusters:
        probas.append(base_pair_proba(cluster, structs, base_pairs_list, fcs, n,s))
        #probas.append(base_pair_proba_v2(cluster, structs, base_pairs_list, n,s))
        plist=[RNA.ep(k[0],k[1],v) for k,v in probas[-1].items()]
        reps.append(RNA.MEA_from_plist(plist,sequences[0], gamma)) #Does the sequence have an influence here?
      
        
    #Clusters' diversity
    clusters_div = [cluster_diversity(i, probas) for i in range(k)]

    #Clusters' ensemble energy
    clusters_ensemble_e = [cluster_ensemble_energy(cluster, structs,  T,fcs, n,s) for cluster in clusters]

    #Clusters' minimum free energy
    clusters_min_e = [min( [fcs[i// n].eval_structure(structs[i]) for i in cluster] ) for cluster in clusters]

    #Printing results

    data = [[len(clusters[i]),clusters_min_e[i],clusters_div[i],clusters_ensemble_e[i],reps[i][0]] for i in range( k)]
    df = pd.DataFrame(data, columns=["Cluster size","Cluster min energy","Cluster diversity","Cluster ensemble energy","MEA representative structure"])
    df.index.name = "Cluster index"
    return df

#Other way:
""" print(tabulate([[i,len(clusters[i]),clusters_min_e[i],clusters_div[i],clusters_ensemble_e[i],reps[i][0]] for i in range( k)], 
        headers=["Cluster index","Cluster size","Cluster min energy","Cluster diversity","Cluster ensemble energy","MEA representative structure"])) """

    
def main(args):
    
    sequences = list(RNA.file_msa_read(args.filename)[2])

    cl_results = clustering(sequences,args.k, args.n)

    df = analyze_clusters(cl_results[0],cl_results[1],cl_results[2],cl_results[3],cl_results[4],args.n,sequences, args.k,args.T, args.gamma)

    print(df)

    if args.outcsv!=None:
        df.to_csv(args.outcsv+".csv")

    if args.outxlsx!=None:
        df.to_excel(args.outxlsx+".xlsx")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename",type=str, help="Name of the Stockholm file representing the RNA alignment (with path)")
    parser.add_argument("-n",type=int, help="Number of generated structures for each sequence of the alignment (default: 10)",default= 30)
    parser.add_argument("-k",type=int, help="Number of clusters (default: 5)",default=5)
    parser.add_argument("-T",type=int, help="Temperature for the computation of the cluster ensemble energy (default: 310.15)",default=310.15)
    parser.add_argument("-gamma",type=int, help="Value of the gamma constant for the computation of the MEA structure (default: 5)",default=5)
    parser.add_argument("-outcsv",type=str, help="Store the resulting DataFrame in a CSV file with the given filename")
    parser.add_argument("-outxlsx",type=str, help="Store the resulting DataFrame in a XLSX file with the given filename")

    args = parser.parse_args()

    main(args)


