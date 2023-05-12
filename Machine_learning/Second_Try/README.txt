We use the same ML script to train and test the model as the first try (exception made on how we make the dummies for the one_hot encoding)

To obtains the data however we used a radically different approch.
We consider that the node we are evaluating (ie. the row we are filling) only has backbone link.
We will now retrieve the datas on his 2 (or 1 for end sequence nucleotide) neighbors (their type, and the number of link of each type they have with their surroundings).
We then proceed to do this for their neighbors (and on and on depending on the DISTANCE value we put at the begining of the script)
We have made sure that the nuclotide we are processing is not being visited again (this would induced a huge biais and not reflect the data we would collect in real life).

Here's a quick explaination of each row of the datas:

reponse: the anwer we want to predict (not kept in the matrice during the ML algorithm for obvious reasons)

at_-1: type of nucleotide that came before him in the sequence ( this node and the "response node" have a B53 edge between them in this specific order).
at_1: type of nucleotide that came after him in the sequence ( the "response node" and this node have a B53 edge between them in this specific order).
=> the values can be ["A","U","G","C","Mod" (for modified nucleotide), "None" (if their is no neighbor) ]
nb_voisin: number of neighbors (1 or 2 depending on where the nculeotide is)

then we have the datas we obtains in a recurssive way:
Given a distance from the response node X we have:
nb_voisin_at_{X}: the number of neighbor at distance X from the successor side
nb_voisin_at_-{X}: the number of neighbor at distance X from the predecessor side

[nuc]_at_{X}: the number of nuclotide of type [nuc] at distance X from the successor side
[nuc]_at_-{X}: the number of nuclotide of type [nuc] at distance X from the predecessor side
nuc can take multiple values:  ["A","U","G","C","Mod"]

[liaison_type]_at_{X}: the number of edge type (link in biological term) that the nucleotide at distance X have from the successor side
[liaison_type]_at_-{X}: the number of edge type (link in biological term) that the nucleotide at distance X have from the predecessor side

[liaison_type] can take multiple values:
['CWW','TWW','CSS','TSS','CHH','THH','CWS','TWS','CWH','TWH','CHS','THS','B53']


FILE DESCRIPTOR:

-get_data_second_try.py
script used to get the node.csv datafile (inside node.zip)
-node2.zip
Dataset used
-subsampling2.py
script used for each different mode (only little modification needed to be made for each optimization

-score_testing2.csv
scoring result of
-grid_cv2.pickle
GridSearchCV (sklearn) object for SVM
-grid_KN2.pickle
GridSearchCV (sklearn) object for KNearestClassification 
-grid_RF2.pickle
GridSearchCV (sklearn) object for RandomForstClassification 











