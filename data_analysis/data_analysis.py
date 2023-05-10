import networkx
import pickle
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import collections as cl
import csv
import pandas as pd


#TO LOAD DATAS

nuc="nucleotide"
regular_list=["A","C","G","U"]
all_nuc=set()

glob=[]
compteur_normal={"A":0,"C":0,"G":0,"U":0}

os.chdir(os.path.join("..","stage_M1","test","graph","NetworkxGraph"))
compteur=0
useless=False
for name in os.listdir():
    if(useless):
        continue
    compteur+=1
    print(compteur,"/6631")
    with open(name, 'rb') as f:
        rna_graph = pickle.load(f)
        f.close()

    info=[]
    for i in rna_graph.nodes:
        print(len(all_nuc))
        if(len(all_nuc)>=326):
            useless=True
            continue
        all_nuc.add((rna_graph.nodes[i][nuc],rna_graph.nodes[i]["nucleotide_one_letter"]))
        if(rna_graph.nodes[i][nuc] not in regular_list):
            import time
            time.sleep(20)
            link=0
            strange=[]
            for voisin in rna_graph.neighbors(i):
                link+=1
                if rna_graph.nodes[voisin][nuc] not in regular_list:
                    strange.append((rna_graph.nodes[voisin][nuc],rna_graph.edges[i,voisin,0]['label']))
            info.append((rna_graph.nodes[i][nuc],link,strange))
        else:
            compteur_normal[rna_graph.nodes[i][nuc]]+=1
    del rna_graph
    glob.append((name,info))

print(len(all_nuc))
jpp={}
for a,b in all_nuc:
    jpp[a]=b

jpp=dict(sorted(jpp.items(), key=lambda item: item[1]))
print(jpp.values())

with open('../parent_nucs.csv', 'w') as f:
    for key in jpp.keys():
        f.write("%s,%s\n"%(key,jpp[key]))


#allows to upload data as pickle
with open(os.path.join("..","info_global.pickle"), 'wb') as handle:
    pickle.dump(glob, handle)
with open(os.path.join("..","list_nuc.pickle"), 'wb') as handle:
    pickle.dump(all_nuc, handle)
with open(os.path.join("..","compteur_normal.pickle"), 'wb') as handle:
    pickle.dump(compteur_normal, handle)


#reloading of data
os.chdir(os.path.join(".","stage_M1","test","graph","NetworkxGraph"))

with open(os.path.join("..","info_global.pickle"), 'rb') as handle:
    glob=pickle.load(handle)
with open(os.path.join("..","list_nuc.pickle"), 'rb') as handle:
    list_nuc=pickle.load(handle)
with open(os.path.join("..","compteur_normal.pickle"), 'rb') as handle:
    compteur_normal=pickle.load(handle)


print(list_nuc)
#MODIFIED NUC CSV

dic={}
for i in list_nuc:
    dic[i]=0
for a,b in glob:
    for nuc,_,_ in b:
        dic[nuc]+=1
dic["A"]=compteur_normal["A"]
dic["C"]=compteur_normal["C"]
dic["G"]=compteur_normal["G"]
dic["U"]=compteur_normal["U"]
print(dic.keys())

dic_sorted = dict(sorted(dic.items(), key=lambda x:x[1]))

with open('../modified_nuc_effectif.csv', 'w') as f:
    for key in dic_sorted.keys():
        f.write("%s,%s\n"%(key,dic_sorted[key]))"


#edge count

linkage=[]
for a,b in glob:
    for nuc,link,_ in b:
        if nuc is not "N":
            linkage.append(link)

print(len(linkage))

compt=cl.Counter(linkage)
with open('../nb_links_wo_N.csv', 'w') as f:
    for key in compt.keys():
        f.write("%s,%s\n"%(key,compt[key]))

print(compt)

#NB edge with specific nucs

linkage=[]
for a,b in glob:
    for nuc,_,link in b:
        if nuc!="N":
            for nuc,type in link:
                if nuc!="N":
                    linkage.append(type)

print(len(linkage))

compt=cl.Counter(linkage)
with open('../nb_links_specific_type_woN.csv', 'w') as f:
    for key in compt.keys():
        f.write("%s,%s\n"%(key,compt[key]))

print(compt)

#detailed edge count
dic_taille={}
for a,b in glob:
    dic_taille[a]=0
    for nuc,_,link in b:
        if nuc is not "N":
            dic_taille[a]+=1


with open('../nb_nuc_per_graph_detailed_woN.csv', 'w') as f:
    for key in dic_taille.keys():
        f.write("%s,%s\n"%(key,dic_taille[key]))


#get the original nucleotide 
nuc_pere={}
with open('./test/graph/parent_nucs.csv', mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            nuc_pere[row['V3']]=row['V1']
nuc_pere['A']='A'
nuc_pere['U']='U'
nuc_pere['C']='C'
nuc_pere['G']='G'
nuc_pere['DG']='DG'
nuc_pere['DU']='DU'
nuc_pere['DC']='DC'
nuc_pere['DA']='DA'


df = pd.read_csv('./test/graph/modified_nuc_effectif.csv', sep=',',header=None)
lig_ids=df[0]
faulty=[]
for a in lig_ids:
    if(a not in nuc_pere):
        faulty.append(a)
print(faulty)