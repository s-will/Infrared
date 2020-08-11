import argparse
import random
import RNA

def remove_gaps(args):
    msa=RNA.file_msa_read(args.infile)
    Id = msa[3]
    

    if args.n==-1:
        aln=msa[2]
        Names = msa[1]
    else:
        indices = random.sample([i for i in range(msa[0])],args.n)
        aln = [msa[2][i] for i in indices]
        Names = [msa[1][i] for i in indices]

    gap_cols = []
    for seq in aln:
        for i in range(len(seq)):
            if seq[i]=='-':
                gap_cols.append(i)

    new_aln=[ ''.join([seq[i] for i in range(len(seq)) if i not in gap_cols]) for seq in aln ]
    
    RNA.file_msa_write(args.outfile,Names,new_aln,Id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Removing gaps from a MSA')

    parser.add_argument('infile', help="Input Stockholm file of the alignment")

    parser.add_argument('outfile', help="Output Stockholm filename of the MSA without gaps")

    parser.add_argument('-n', type=int,default=-1,help="Size of the alignment to return (subset of initial one)")
    
    args=parser.parse_args()

    remove_gaps(args)