from Bio import SeqIO
import sys

products=["DNA gyrase subunit A", "DNA gyrase subunit B", "DNA topoisomerase (ATP-hydrolyzing) subunit B", "DNA topoisomerase IV subunit B", "DNA topoisomerase IV subunit A"]

genes = ["gyrA", "gyrB", "parE", "parC"]

corres={"DNA gyrase subunit A":"gyrA", "DNA gyrase subunit B":"gyrB", "DNA topoisomerase (ATP-hydrolyzing) subunit B": "gyrB", "DNA topoisomerase IV subunit B":"parE", "DNA topoisomerase IV subunit A":"parC"}

def extract_translation(filename, gen, prod, corr):
    strain = ""
    for record in SeqIO.parse(filename, "genbank"):
        for f in record.features:
            if f.type=="source":
                try:
                    strain = f.qualifiers["strain"][0].replace(" ", "_")
                except KeyError:
                    try:
                        strain = f.qualifiers["isolate"][0].replace(" ", "_")
                    except KeyError:
                        print(f.qualifiers)
            if f.type=="CDS":
                try:
                    if f.qualifiers["gene"][0] in gen:
                        with open(f.qualifiers["gene"][0]+".fasta", "a") as myfile:
                            try:
                                header = ">"+record.annotations["organism"].replace(" ", "_")+"_"+strain+"_"+f.qualifiers["locus_tag"][0]+"\n"
                                myfile.write(header)
                                myfile.write(f.qualifiers["translation"][0]+"\n")
                            except KeyError:
                                print(f)
                except KeyError:
                    if f.qualifiers["product"][0] in prod:
                        with open(corr[f.qualifiers["product"][0]]+".fasta", "a") as myfile:
                            try:
                                header = ">"+record.annotations["organism"].replace(" ", "_")+"_"+strain+"_"+f.qualifiers["locus_tag"][0]+"\n"
                                myfile.write(header)
                                myfile.write(f.qualifiers["translation"][0]+"\n")
                            except KeyError:
                                print(f)
    
extract_translation(sys.argv[1], genes, products, corres) 
