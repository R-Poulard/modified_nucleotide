import pickle
import csv, sys
import os
import os.path
from pandas import *
import seaborn as sns
import pandas as pd
import rdkit.Chem as Chem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
import numpy
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#Script to retrieve SMILES comes from here and have been modified
#https://github.com/MunibaFaiza/cheminformatics/blob/main/pdb_ligand_id-to-smi.ipynb


# Reading ligand IDs
df = pd.read_csv('./test/graph/modified_nuc_effectif.csv', sep=',',header=None)

# defining output dir
lig_dir = "./test/ligands-sdf"



#os.mkdir(lig_dir)
os.chdir(lig_dir)

lig_ids=df[0]
print(lig_ids)
print(len(lig_ids))


for id in lig_ids:
   cmd = os.system("wget "+" "+"https://files.rcsb.org/ligands/download/"+id+"_ideal.sdf")

join_cmd = os.system("cat *.sdf > ligands.sdf")
joined_ligands = "ligands.sdf"

def converter(file_name):
    
    sdf_file = Chem.SDMolSupplier(file_name)
    listing=[]
    compteur=0
    for mol in sdf_file:
        if mol is not None:             # avoiding compounds that cannot be loaded.
            smi = Chem.MolToSmiles(mol)
            listing.append((lig_ids[compteur],f"{smi}"))
        else:
            print("MISSING DATA",lig_ids[compteur])
        compteur+=1
    return listing


if __name__ == "__main__":
    l=converter(joined_ligands)
    #Manually add 6FU beceause MolToSmiles seems to not work for him (don't know why)
    mol="c1cn(c(=O)nc1N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(O)(O)O)O)SC(F)(F)F"
    smi = Chem.MolFromSmiles(mol)
    l.append(('6FC',Chem.MolToSmiles(smi)))

    #creation of smiles dictionnay and maccs dictionnary
    dic={}
    for a,b in l:
        dic[a]=b
    with open(os.path.join("..","dic_smiles.pickle"), 'wb') as handle:
        pickle.dump(dic, handle)
    dic_maccs={}
    for a,b in dic.items():
        dic_maccs[a]= MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(b))
    with open(os.path.join("..","dic_maccs_unfolded.pickle"), 'wb') as handle:
        pickle.dump(dic_maccs, handle)

    import csv

    #HEATMAP CREATION
    #reading original nucleotide to order the similarity heatmap by parents
    nuc_pere={}
    with open('../graph/parent_nucs.csv', mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            nuc_pere[row['V3']]=row['V1']
    
    """
    #USED FOR CLOSER HEATMAP
    dic_maccs={"A":dic_maccs["A"],"G":dic_maccs["G"],"C":dic_maccs["C"],"U":dic_maccs["U"],"PSU":dic_maccs["PSU"],"5MC":dic_maccs["5MC"],"OMG":dic_maccs["OMG"],"5MU":dic_maccs["5MU"],
    "2MG":dic_maccs["2MG"],"OMC":dic_maccs["OMC"],"A2M":dic_maccs["A2M"]
    }
    """
    mapping=numpy.array([[0.0]*len(dic_maccs)]*len(dic_maccs))
    i_ind=0
    for i in dic_maccs:
        print(i_ind)
        y_ind=0
        for y in dic_maccs:
            mapping[i_ind,y_ind]=DataStructs.FingerprintSimilarity(dic_maccs[i], dic_maccs[y])
            y_ind+=1
        i_ind+=1
    print(mapping)
    df=pd.DataFrame(mapping,columns=list(dic_maccs.keys()),index=list(dic_maccs.keys()))
    df.to_csv("finger_print_corr.csv")
    sns.set(font_scale=1)
    ht=sns.heatmap(df)
    plt.show()

    #Create dendogram from correlation 
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import squareform

    plt.figure(figsize=(12,5))
    dissimilarity = 1 - abs(df)
    Z = linkage(squareform(dissimilarity), 'complete')

    dendrogram(Z, orientation='top', 
           leaf_rotation=90);
    plt.show()

    sns.clustermap(df)
    plt.show()
