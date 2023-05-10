import networkx
import pickle
import os
import matplotlib
import matplotlib.pyplot as plt
import collections as cl
import csv
import pandas as pd
"""
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
        f.write("%s,%s\n"%(key,jpp[key]))"""
"""

with open(os.path.join("..","info_global.pickle"), 'wb') as handle:
    pickle.dump(glob, handle)
with open(os.path.join("..","list_nuc.pickle"), 'wb') as handle:
    pickle.dump(all_nuc, handle)
with open(os.path.join("..","compteur_normal.pickle"), 'wb') as handle:
    pickle.dump(compteur_normal, handle)"""

"""
os.chdir(os.path.join(".","stage_M1","test","graph","NetworkxGraph"))

with open(os.path.join("..","info_global.pickle"), 'rb') as handle:
    glob=pickle.load(handle)
with open(os.path.join("..","list_nuc.pickle"), 'rb') as handle:
    list_nuc=pickle.load(handle)
with open(os.path.join("..","compteur_normal.pickle"), 'rb') as handle:
    compteur_normal=pickle.load(handle)
"""
"""
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
print(dic.keys())"""

"""
dic_sorted = dict(sorted(dic.items(), key=lambda x:x[1]))


with open('../modified_nuc_effectif.csv', 'w') as f:
    for key in dic_sorted.keys():
        f.write("%s,%s\n"%(key,dic_sorted[key]))"""


#linkage compteur

"""linkage=[]
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
"""

"""
#NB LINK SPECIFIQUE

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

"""

"""#NB LINK SPECIFIQUE

linkage=[]
for a,b in glob:
    for nuc,_,link in b:
        for nuc,type in link:
            if nuc!="N":
                linkage.append(type)

print(len(linkage))

compt=cl.Counter(linkage)
with open('../nb_links_specific_type_woN.csv', 'w') as f:
    for key in compt.keys():
        f.write("%s,%s\n"%(key,compt[key]))

print(compt)


dic_taille={}
for a,b in glob:
    dic_taille[a]=0
    for nuc,_,link in b:
        if nuc is not "N":
            dic_taille[a]+=1


with open('../nb_nuc_per_graph_detailed_woN.csv', 'w') as f:
    for key in dic_taille.keys():
        f.write("%s,%s\n"%(key,dic_taille[key]))"""

nuc_pere={}
with open('./parent_nucs.csv', mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            nuc_pere[row['V3']]=row['V1']
nuc_pere['A']='A'
nuc_pere['U']='U'
nuc_pere['C']='C'
nuc_pere['G']='G'
nuc_pere['DG']='G'
nuc_pere['DU']='U'
nuc_pere['DC']='C'
nuc_pere['DA']='A'


df = pd.read_csv('./modified_nuc_effectif.csv', sep=',',header=None)
lig_ids=df[0]
faulty=[]
for a in lig_ids:
    if(a not in nuc_pere):
        faulty.append(a)

nuc="nucleotide"
regular_list=["A","C","G","U","N"]
link_list=['CWW','TWW','CSS','TSS','CHH','THH','CWS','TWS','CWH','TWH','CHS','THS','B53']


os.chdir(os.path.join(os.sep,"srv","rna","NetworkxGraph"))
compteur=0
full_count=0
info={}
info['reponse']=[]
info['parent']=[]
info['nb_voisin']=[]
info['nb_voisin_mod']=[]

for lab in link_list:
    info[lab]=[]
for name in os.listdir():
    
    compteur+=1
    print(compteur,"/6631")
    with open(name, 'rb') as f:
        rna_graph = pickle.load(f)
        f.close()

    for i in rna_graph.nodes:
        full_count+=1
        if rna_graph.nodes[i][nuc] not in faulty:
            
            link=0
            special_link=0
            tmp={}
            for lab in link_list:
                tmp[lab]=0
            for voisin in rna_graph.neighbors(i):
                link+=1
                if rna_graph.nodes[voisin][nuc] not in regular_list:
                    special_link+=1
                for a,b in rna_graph.get_edge_data(i,voisin).items():
                    
                    if b['label'] in link_list and b['near']==False:
                        tmp[b['label']]+=1
            info['reponse'].append(rna_graph.nodes[i][nuc])   
            info['nb_voisin'].append(link)
            info['nb_voisin_mod'].append(special_link)
            info['parent'].append(nuc_pere[rna_graph.nodes[i][nuc]])
            for i in link_list:
                info[i].append(tmp[i])
    del rna_graph

for i in info:
    print(i,len(info[i]))
data=pd.DataFrame(info)
print(data.shape)
data.to_csv(os.path.join("~","M1","get_nodes","node.csv"),header=True,index=False)
