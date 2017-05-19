#!/usr/bin/env python
import Bio
import Bio.SeqIO
import Bio.AlignIO
import sys

def read_alignment(filename):
    return(Bio.AlignIO.parse(filename, "fasta"))


def return_seqRec(aliRec):
    seqs={}
    for i in aliRec:
        for j in i:
            seqs[j.id]=j.seq
    return(seqs)

def return_avg_identity(seqRec):
    print("\t"+"\t".join(seqRec.keys()))
    for i in seqRec.keys():
        id=i+"\t"
        for j in seqRec.keys():
            matches=sum(aa1 == aa2 and aa1 != "-" for aa1, aa2 in zip(seqRec[i], seqRec[j]))
            unmatches=sum(aa1 == "-" for aa1 in seqRec[i])
            id+=str(matches/float(len(seqRec[i])-unmatches))+" ("+str(len(seqRec[i])-unmatches)+")\t"
        print(id)
        
        


if __name__ == "__main__":
    print(sys.argv[1])
    return_avg_identity(return_seqRec(read_alignment(sys.argv[1])))
