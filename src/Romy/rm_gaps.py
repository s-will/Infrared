#!/usr/bin/env python3

import argparse
import random
import RNA
import rna_support 

def remove_gaps(args):
    msa=RNA.file_msa_read(args.infile)
    Id = msa[3]
    consensus_struct=msa[4]

    if args.n==-1:
        n = msa[0]
        aln=msa[2]
        Names = msa[1]
    else:
        n = args.n
        indices = random.sample([i for i in range(msa[0])],args.n)
        aln = [msa[2][i] for i in indices]
        Names = [msa[1][i] for i in indices]
        
    
    N = len(aln[0])

    gap_cols = [False for _ in range(N)]
    for i in range(N):
        for seq in aln:
            if seq[i]=='-':
                gap_cols[i]=True
                break


    new_aln=[ ''.join([seq[i] for i in range(N) if not gap_cols[i]]) for seq in aln ]
    
    pairing=rna_support.parseRNAStructure(consensus_struct)
    for x in [i for i in range(N) if gap_cols[i]]:
        if pairing[x]!=-1:
            consensus_struct=consensus_struct[:pairing[x]]+'.'+consensus_struct[pairing[x]+1:]    
    new_struct=''.join([consensus_struct[i] for i in range(N) if not gap_cols[i]])

    
    if args.outfile!=None:
        RNA.file_msa_write(args.outfile,Names,new_aln,Id,new_struct,'',2)

    return [n,Names,aln]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Removing gaps from a MSA')

    parser.add_argument('infile', help="Input Stockholm file of the alignment")

    parser.add_argument('-outfile',help="Output Stockholm filename of the MSA without gaps")

    parser.add_argument('-n', type=int,default=-1,help="Size of the alignment to return (subset of initial one)")
    
    args=parser.parse_args()

    remove_gaps(args)