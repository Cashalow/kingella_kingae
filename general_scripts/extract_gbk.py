from Bio import SeqIO
import sys

products=["DNA gyrase subunit A", "DNA gyrase subunit B", "DNA topoisomerase (ATP-hydrolyzing) subunit B", "DNA topoisomerase IV subunit B"]
genes = ["gyrA", "gyrB", "parE"]

corres={"DNA gyrase subunit A":"gyrA", "DNA gyrase subunit B":"gyrB", "DNA topoisomerase (ATP-hydrolyzing) subunit B": "gyrB", "DNA topoisomerase IV subunit B":"parE"}

def extract_translation(filename, gen, prod, corr):
    for record in SeqIO.parse(filename, "genbank"):
        for f in record.features:
            if f.type=="CDS":
#                print(f.qualifiers["product"])
                try:
                    if f.qualifiers["gene"][0] in gen:
                        with open(f.qualifiers["gene"][0]+".fasta", "a") as myfile:
                            try:
                                myfile.write(">"+record.annotations["organism"].replace(" ", "_")+"_"+f.qualifiers["locus_tag"][0]+"\n")
                                myfile.write(f.qualifiers["translation"][0]+"\n")
                            except KeyError:
                                print(f)
                except KeyError:
                    if f.qualifiers["product"][0] in prod:
                        with open(corr[f.qualifiers["product"][0]]+".fasta", "a") as myfile:
                            try:
                                myfile.write(">"+record.annotations["organism"].replace(" ", "_")+"_"+f.qualifiers["locus_tag"][0]+"\n")
                                myfile.write(f.qualifiers["translation"][0]+"\n")
                            except KeyError:
                                print(f)
    
extract_translation(sys.argv[1], genes, products, corres) 
