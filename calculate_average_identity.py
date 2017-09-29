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
#            print(i)
 #           print(j)
            nogaps = sum(aa1 != "-" and aa2 != "-" for aa1, aa2 in zip(seqRec[i], seqRec[j]))
            matches = sum(aa1 == aa2 and aa1 != "-" and aa1 != "n" and aa1 != "N" for aa1, aa2 in zip(seqRec[i], seqRec[j]))
            seqlen = sum(aa1 != "-" for aa1 in seqRec[i])
            aln = sum(aa1 != "-" and aa2 != "-" for aa1, aa2 in zip(seqRec[i], seqRec[j]))
            id+="{:2.2f} ({:2.2f}%)\t".format(matches/float(nogaps)*100, aln/float(seqlen)*100)
        print(id)
        
        


if __name__ == "__main__":
    print(sys.argv[1])
    return_avg_identity(return_seqRec(read_alignment(sys.argv[1])))
