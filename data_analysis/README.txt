Data analysis, couldn't put the data on here beceause the archive is 3,1Go

FILE DESCRIPTOR:

-atom_mappings_refined.txt
text file with modified_nucleotide[tabspace]original nuceotide 

-deta_analysis.py
python script used to obtain all those files

-detail_graph_woN.png
histogramm of the number of graph with more that 1 modified nuceotide 

-graph_detailed_woN_sub50.png
histogram of the number of graphi wtih less than 50 modified nucleotide

-info_global.pickle
pickle file use to obtain all the infos 
list of tuple :
(graph name, information)
with information as a list of tuple:
(name of modified nucleotide, number of neighbors, list of (label edge, modified nucleotid neighbors)
*sorry for this format*

-list_nuc.pickle
set of all the nuclotide type

-nb_links_specific_type_woN.csv
number of edge each nucleotide have with other modified nucleotid (gathered on a table)

-nb_links_woN.png
barplot of the number of edge modified nucleotid have (segregation done on edge label)

-nb_nuc_per_graph_detailed_woN.csv
csv file with each row = name_graph,number of modified nucleotide

-parent_nucs.csv
atom_mappings_refined.txt as csv
