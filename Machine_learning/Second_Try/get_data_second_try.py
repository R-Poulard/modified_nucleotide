import networkx
import pickle
import os
import matplotlib
import matplotlib.pyplot as plt
import collections as cl
import csv
import pandas as pd

DISTANCE=1
link_list=['CWW','TWW','CSS','TSS','CHH','THH','CWS','TWS','CWH','TWH','CHS','THS','B53']
regular_list=["A","C","G","U"]
nuc="nucleotide"


def create_main_dic(distance):
    tmp={}
    tmp['reponse']=[]
    tmp['at_-1']=[]
    tmp['at_1']=[]
    tmp['nb_voisin']=[]
    iter= range(1,distance+1)
    for index in iter:
        tmp['nb_voisin_at_'+str(index)]=[]
        for value in ["Mod","A","C","G","U"]:
            tmp[value+"_at_"+str(index)]=[]
        for edg in link_list:
            tmp[edg+"_at_"+str(index)]=[]
    iter=range((distance)*-1,0,1)
    for index in iter:
        tmp['nb_voisin_at_'+str(index)]=[]
        for value in ["Mod","A","C","G","U"]:
            tmp[value+"_at_"+str(index)]=[]
        for edg in link_list:
            tmp[edg+"_at_"+str(index)]=[]
    return tmp

def create_fake_list(distance,index,successor):
    tmp=[]
    if successor:
         iter=range(index+1,distance+1)
    else:
         iter=range(-distance,index,1)
    for index in iter:
         tmp.append('nb_voisin_at_'+str(index))
         for value in ["Mod","A","C","G","U"]:
            tmp.append(value+"_at_"+str(index))
         for edg in link_list:
            tmp.append(edg+"_at_"+str(index))
    return tmp
def create_fake(distance,successor=True):
     tmp={}
     if successor:
          iter= range(1,distance+1)
     else:
          iter=range((distance)*-1,0,1)
     for index in iter:
             tmp['nb_voisin_at_'+str(index)]=0
             for value in ["Mod","A","C","G","U"]:
                 tmp[value+"_at_"+str(index)]=0
             for edg in link_list:
                 tmp[edg+"_at_"+str(index)]=0 
     return tmp


def get_surrounding(G,not_count,source,successor=True,last=False,index=0):
      tmp={}
      need=True
      tmp['nb_voisin_at_'+str(index)]=0
      for value in ["Mod","A","C","G","U"]:
         tmp[value+"_at_"+str(index)]=0
      for edg in link_list:
         tmp[edg+"_at_"+str(index)]=0
      if successor:
         iter=G.successors(source)
      else:
         iter=G.predecessors(source)
      for voisin in iter:
           if voisin not in not_count:
             try:
                tmp[G.nodes[voisin][nuc]+"_at_"+str(index)]+=1
             except KeyError:
                tmp["Mod_at_"+str(index)]+=1
             if successor:
                 dic_edge=rna_graph.get_edge_data(source,voisin)
             else:
                 dic_edge=rna_graph.get_edge_data(voisin,source)
             for a,b in dic_edge.items():
                if b['near']==False and b['label'] in link_list:
                   tmp['nb_voisin_at_'+str(index)]+=1
                   tmp[b['label']+"_at_"+str(index)]+=1
                   if not last:
                       need=False
                       if successor:
                           surr=get_surrounding(G,not_count+[source],voisin,successor,last=(index+1==DISTANCE),index=index+1)
                       else:
                           #print("pile",index,-DISTANCE)
                           surr=get_surrounding(G,not_count+[source],voisin,successor,last=(index-1==-DISTANCE),index=index-1)
                       for key in surr:
                          if key not in tmp:
                             tmp[key]=surr[key]
                          else:
                             tmp[key]+=surr[key]
      if need:
         fake=create_fake_list(DISTANCE,index,successor)
         for key in fake:
             tmp[key]=0
      return tmp
"""
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
"""

os.chdir(os.path.join(".","stage_M1","test","ML","graph_test"))
compteur=0

info=create_main_dic(DISTANCE)

for name in os.listdir():
    compteur+=1
    print(compteur,"/6631")
    with open(name, 'rb') as f:
        rna_graph = pickle.load(f)
        f.close()
    for i in rna_graph.nodes:            
            index_at_one=-1
            index_at_zero=-1
            link=0
            tmp={}
            tmp["at_-1"]="None"
            tmp["at_1"]="None"
            for lab in link_list:
                tmp[lab+"_at_-1"]=0
                tmp[lab+"_at_1"]=0
            for lab in regular_list:
                tmp[lab+"_at_-1"]=0
                tmp[lab+"_at_1"]=0
            print("pour ",i)
            for voisin in rna_graph.successors(i):
                    print(voisin)
                    for a,b in rna_graph.get_edge_data(i,voisin).items():
                        if b['label']=='B53':
                            link+=1
                            index_at_one=voisin
                            if rna_graph.nodes[voisin][nuc] in regular_list:
                                tmp["at_1"]=rna_graph.nodes[voisin][nuc]
                            else:
                                tmp["at_1"]="Mod"
            for voisin in rna_graph.predecessors(i):
                    for a,b in rna_graph.get_edge_data(voisin,i).items():
                       if b['label']=='B53':
                            link+=1
                            index_at_zero=voisin
                            if rna_graph.nodes[voisin][nuc] in regular_list:
                           	    tmp["at_-1"]=rna_graph.nodes[voisin][nuc]
                            else:
                                tmp["at_-1"]="Mod"
            tmp['nb_voisin']=link
            tmp['reponse']=rna_graph.nodes[i][nuc]
            if tmp["at_-1"]!="None":
                dic_before=get_surrounding(rna_graph,[i],index_at_zero,successor=False,last=(-1==-DISTANCE),index=-1)
            else:
                dic_before=create_fake(DISTANCE,successor=False)
            if tmp["at_1"]!="None":
                dic_after=get_surrounding(rna_graph,[i],index_at_one,successor=True,last=(1==DISTANCE),index=1)
            else:
                dic_after=create_fake(DISTANCE,successor=True)
            for key in dic_before:
                #print(key)
                tmp[key]=dic_before[key]
            for key in dic_after:
                #print(key)
                tmp[key]=dic_after[key]
            for key in info:
                info[key].append(tmp[key])
for key in info:
	print(key,len(info[key]),end="\t")
data=pd.DataFrame(info)
print(data)
data.to_csv(os.path.join("..","TRY.csv"),header=True,index=False)
#print("total of suppressed edges",compteur_suppr_edge)
