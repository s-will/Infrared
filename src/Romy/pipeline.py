import romy
import clustering as cl
import rm_gaps as rg
import random
import seaborn as sns
import matplotlib.pyplot as plt 
from scipy import stats


import argparse
import pandas as pd
import RNA

## What this script does:
# 1) Remove gaps of the alignment contained in the infile argument, 
#    and takes a random subset of size -ss 
# 2) Sample -n alignments targeting features of the input alignment.
# 3) For each sampled alignment, we sample structures and cluster them 
#    to store the minimum free energy and ensemble energy of the 2 best clusters.

## The resulting dataframe contains all these results and its last line is
#  for the input alignment.

def main(args):
    args.infile= max_dist_seqs(args.ss,args.infile)
    alignments=romy.main(args)
    if args.ss == -1: ss = len(alignments[0])
    else: ss = args.ss
    x = int(1000/ss) #Number of structures to generate per sequence to have 1000 structures for the alignment
    res= {"best_min_e":[], "best_ens_e" :[], "second_best_min_e" :[], "second_best_ens_e" :[]}
    for sequences in alignments:
        #print(sequences)
        cl_results = cl.clustering(sequences,args.beta,args.k,x)
        df = cl.analyze_clusters(cl_results[0],cl_results[1],cl_results[2],cl_results[3],cl_results[4],cl_results[5],cl_results[6],x,sequences,args.k,args.T,args.gamma)
        i1, i2 = minimin(df)
        res["best_min_e"].append(df["Cluster min energy"][i1])
        res["second_best_min_e"].append(df["Cluster min energy"][i2])
        res["best_ens_e"].append(df["Cluster ensemble energy"][i1])
        res["second_best_ens_e"].append(df["Cluster ensemble energy"][i2])
    
    #Then collecting results from the initial alignments (gapless)
    sequences = list(RNA.file_msa_read(args.infile)[2])
    msa_res = cl.clustering(sequences,args.beta,args.k,x)
    df = cl.analyze_clusters(msa_res[0],msa_res[1],msa_res[2],msa_res[3],msa_res[4],msa_res[5],msa_res[6],x,sequences,args.k,args.T,args.gamma)
    i1, i2 = minimin(df)
    res["best_min_e"].append(df["Cluster min energy"][i1])
    res["second_best_min_e"].append(df["Cluster min energy"][i2])
    res["best_ens_e"].append(df["Cluster ensemble energy"][i1])
    res["second_best_ens_e"].append(df["Cluster ensemble energy"][i2])
    
    
    results = pd.DataFrame.from_dict(res)
    results.to_csv(args.outcsv+".csv")

    if args.plot_res:
        #Compare the distribution 
        dist = results["second_best_min_e"][:-1]
        point = [ results["second_best_min_e"][len(results)-1] ]
        print("Statistics: ",stats.ks_2samp(dist, point))
        
        plt.plot([results["second_best_min_e"][len(results)-1]],[0],marker='o',markersize = 10, color='red')
        sns.distplot(results["second_best_min_e"][:-1],hist=True, kde=True, 
             color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})
        

        plt.title("Sampled distribution of the 2nd best cluster min energy")
        plt.show()
    results.to_csv(args.outcsv+".csv")

    return results

def minimin(df):
    m1, m2 = float("inf"),float("inf")
    i1, i2 = -1, -1
    for i in range(len(df)):
        x = df["Cluster min energy"][i]
        if x<m1:
            m1,m2 = x,m1
            i1,i2 = i, i1
        elif x<m2:
            m2 = x
            i2 = i
    
    return i1,i2

def max_dist_seqs(subset_size,infile):

    msa=RNA.file_msa_read(args.infile)
    sequences = list(msa[2])
    if subset_size==-1 or subset_size>len(sequences):
        return infile
    else:
        Id = msa[3]
        consensus_struct=msa[4]
 
        #Compute matrix of hamming distances between sequences
        #Cluster them to have subset_size clusters and choose one rep. from each
        n = len(sequences)
        dist_matrix= [[0 for _ in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i,n):
                v = romy.hamming_distance(sequences[i],sequences[j])
                dist_matrix[i][j],dist_matrix[j][i]=v,v
    
    kmeans = cl.cluster_kmeans(subset_size,dist_matrix)
    clusters = [[i for i in range(n) if kmeans.labels_[i]==k] for k in range(subset_size)]
    indices = [random.choice(clust) for clust in clusters]
    new_aln = [sequences[i] for i in indices]
    names = [msa[1][i] for i in indices]
    outfile = infile.split('.')[0]+'-max-dist-subset.'+infile.split('.')[1]
    RNA.file_msa_write(outfile,names,new_aln,Id,consensus_struct,'',2)
    
    return outfile
    

                

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Boltzmann sampling of homologous sequences')
    parser.add_argument('infile',type=str, help="Input Stockholm file of the alignment")
    parser.add_argument('--ss', type=int, default=-1, help="Subset size: size of the alignment subset to do the analysis on")
    parser.add_argument('--plot_res',action="store_true",help="Plot the resulting sampled distribution and the true value")
    #Arguments for romy sampler
    parser.add_argument('--n','--number', type=int, default=100, help="Number of samples")
    parser.add_argument('--td', type=str, default="nx",help="Method for tree decomposition (see --listtds)")
    parser.add_argument('--listtds', action="store_true", help="List available tree decomposition methods")
    parser.add_argument('--seed', type=int, default=None,
                        help="Seed infrared's random number generator (def=auto)")
    parser.add_argument("--newick",type=str,default=None,
                        help="Filename of the newick phylogenetic tree to use")       
    parser.add_argument('--struct', type=str, default=None, help="Consensus structure for the alignment")
    parser.add_argument('--gc_tolerance', type=float, default=10, help="Target tolerance for the GC content")
    parser.add_argument('--energy_tolerance', type=float, default=10, help="Target tolerance for energies")
    parser.add_argument('--distance_tolerance', type=float, default=4, help="Target tolerance for hamming distances")
    parser.add_argument('--gc_weight', type=float, default=1, help="GC weight")
    parser.add_argument('--energy_weight', type=float, default=1, help="Energy weight")
    parser.add_argument('--distance_weight', type=float, default=1, help="Distance weight")
    parser.add_argument('--plot_td', action="store_true",
                        help="Plot tree decomposition")
    parser.add_argument('--mdbs', action="store_true",
                        help="Perform multi-dim Boltzmann sampling to aim at targets")
    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")
    
    #Arguments for clustering
    parser.add_argument("--k",type=int, help="Number of clusters (default: 5)",default=5)
    parser.add_argument("--T",type=int, help="Temperature for the computation of the cluster ensemble energy (default: 310.15)",default=310.15)
    parser.add_argument("--gamma",type=int, help="Value of the gamma constant for the computation of the MEA structure (default: 5)",default=5)
    parser.add_argument("--beta",type=float, help="BetaScale for the sampling of structures",default= 2.0)
    parser.add_argument("--outcsv",type=str,default="background", help="Store the resulting DataFrame in a CSV file with the given filename")

    args=parser.parse_args()
    main(args)