from bs4 import BeautifulSoup as bs
import csv
import sys
import os
from collections import OrderedDict



def extract_data(filename):
    dic={}
    soup=bs(open(filename).read(), "html.parser")

    for lol in soup.find_all("h5"):
        if lol.get_text() =="Genomes":
            table=lol.find_next("table")
            values=table.find_next("tbody").findAll("tr")
            header=table.find_next("thead").findAll("th")
            l=[k.get_text() for k in header]
            for u in values:
                cells=u.findAll("td")
                v=[s.strip() for s in [k.get_text() for k in cells]]
                dic[v[0]]={}
                for i in range(1,8):
                    dic[v[0]][l[i]]=v[i]

    if not len(dic):
        print(filename)
        raise ValueError("No genome found in this entry")


    for lol in soup.find_all("h5"):
        if lol.get_text() =="Specimen":
            table=lol.find_next("table")
            header=table.find_next("thead").findAll("th")
            l=[k.get_text() for k in header]
            values=table.find_next("tbody").findAll("tr")
            for u in values:
                cells=u.findAll("td")
                v=[s.strip() for s in [k.get_text() for k in cells]]
                if v[0] in dic.keys():
                    dic[v[0]][l[2]]=v[2]



    for lol in soup.find_all("h5"):
        if lol.get_text() =="Drug susceptibility testing":
            table=lol.find_next("table")
            header=table.find_next("thead").findAll("th")
            l=[k.get_text() for k in header]
            values=table.find_next("tbody").findAll("tr")
            for u in values:
                cells=u.findAll("td")
                v=[s.strip() for s in [k.get_text() for k in cells]]
                if v[0] in dic.keys():
                    dic[v[0]][v[1]]={}
                    for i in range(2, len(v)):
                        dic[v[0]][v[1]][l[i]]=v[i]
            
    return(dic)

final_dic={}
files = [f for f in os.scandir("..") if os.path.isfile(f)]
print(files)
for i in files:
    final_dic={** final_dic, ** extract_data("../"+i.name)}
    
     

tests=["Lowenstein-Jensen", "Bactec", "GeneXpert", "Hain", "DST"]


print(final_dic.keys())
for i in final_dic.keys():
    with open('samples.csv', 'a') as f:
        tmp = {}
        tmp["Name"]=i
        tmp = {** tmp, ** {k: final_dic[i][k] for k in list(final_dic[i].keys())[1:8]}}        
        w = csv.DictWriter(f, tmp)
        if os.stat('samples.csv').st_size==0:
            w.writeheader()
        w.writerow(tmp)
        tests_made = [val for val in tests if val in list(final_dic[i].keys())]
        for test in tests_made:
            with open(test+".csv", "a") as g:
                tmpres = {}
                tmpres["Name"]=i
                tmpres = {** tmpres, ** final_dic[i][test]}
                w = csv.DictWriter(g, tmpres)
                if os.stat(test+".csv").st_size==0:
                    w.writeheader()
                w.writerow(tmpres)
        
        
        
