import romy
import clustering as cl
import rm_gaps as rg

import argparse
import pandas as pd


def main(args):
    args.outfile=args.infile.split('.')[0]+'-gapless.'+args.infile.split('.')[1]
    rg.remove_gaps(args)
    args.infile= args.outfile
    alignments,seqnum=romy.main(args)
    x = int(1000/seqnum) #Number of structures to generate per sequence to have 1000 structures for the alignment
    res= {"best_min_e":[], "best_ens_e" :[], "second_best_min_e" :[], "second_best_ens_e" :[]}
    for sequences in alignments:
        cl_results = cl.clustering(sequences,args.k,x)
        df = cl.analyze_clusters(cl_results[0],cl_results[1],cl_results[2],cl_results[3],cl_results[4],x,sequences,args.k,args.T,args.gamma)
        i1, i2 = minimin(df)
        res["best_min_e"].append(df["Cluster min energy"][i1])
        res["second_best_min_e"].append(df["Cluster min energy"][i2])
        res["best_ens_e"].append(df["Cluster ensemble energy"][i1])
        res["second_best_ens_e"].append(df["Cluster ensemble energy"][i2])
    results = pd.DataFrame.from_dict(res)
    print(results)
    if args.outcsv!=None:
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Boltzmann sampling of homologous sequences')
    parser.add_argument('infile',type=str, help="Input Stockholm file of the alignment")
    parser.add_argument('-ss', type=int, default=-1, help="Subset size: size of the alignment subset to do the analysis on")
    
    #Arguments for romy sampler
    parser.add_argument('-n','--number', type=int, default=100, help="Number of samples")
    parser.add_argument('--td', type=str, default="nx",help="Method for tree decomposition (see --listtds)")
    parser.add_argument('--listtds', action="store_true", help="List available tree decomposition methods")
    parser.add_argument('--seed', type=int, default=None,
                        help="Seed infrared's random number generator (def=auto)")
    parser.add_argument("--newick",type=str,default=None,
                        help="Filename of the newick phylogenetic tree to use")       
    parser.add_argument('--struct', type=str, default=None, help="Consensus structure for the alignment")
    parser.add_argument('--gc_tolerance', type=float, default=5, help="Target tolerance for the GC content")
    parser.add_argument('--energy_tolerance', type=float, default=5, help="Target tolerance for energies")
    parser.add_argument('--distance_tolerance', type=float, default=1, help="Target tolerance for hamming distances")
    parser.add_argument('--gc_weight', type=float, default=1, help="GC weight")
    parser.add_argument('--energy_weight', type=float, default=1, help="Energy weight")
    parser.add_argument('--distance_weight', type=float, default=1, help="Distance weight")
    parser.add_argument('--plot_td', action="store_true",
                        help="Plot tree decomposition")
    parser.add_argument('--mdbs', action="store_true",
                        help="Perform multi-dim Boltzmann sampling to aim at targets")
    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")
    
    #Arguments for clustering
    parser.add_argument("-k",type=int, help="Number of clusters (default: 5)",default=5)
    parser.add_argument("-T",type=int, help="Temperature for the computation of the cluster ensemble energy (default: 310.15)",default=310.15)
    parser.add_argument("-gamma",type=int, help="Value of the gamma constant for the computation of the MEA structure (default: 5)",default=5)

    parser.add_argument("-outcsv",type=str, help="Store the resulting DataFrame in a CSV file with the given filename")

    args=parser.parse_args()
    main(args)