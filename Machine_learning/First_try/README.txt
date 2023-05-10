In this first try we present each node with the following attributes
[reponse] what we want to predict
[parent] originel nodes (suppressed for some analysis beceause we show it induced a biais that lead to "too good prediction")
ENCODED IN ONE_HOT DURING THE PROCESSING
[nb_voisins] number of neighbors
[nb_voisins_modifié] number of modified neighbors
then follow 13 attributes on the different edges type the nucleotide have

To predict we use 3 different model (using scklearn)
-SVC
-KNearestClassification
-RandomTreeForst

We use a SearchGridCV each time to find the best parameters (base on f1 criterion) then score using 
-f1_score
-balanced accuracy score
-average precision score
-averege recall score


We tried 3 different optimizaiton:
[1] sub-sampling
[2] no subsampling #not done yet still running (or I gave up as the dataset is really big)
[3] sub-sampling wo the [parent] node

FILE DESCRIPTOR:

-noce.csv
Dataset used
-subsampling.py
script used for each different mode (only little modification needed to be made for each optimization

-score_testing.csv
scoring result of [1]
-grid_cv.pickle
GridSearchCV (sklearn) object for SVM [1]
-grid_KN.pickle
GridSearchCV (sklearn) object for KNearestClassification [1]
-grid_RF.pickle
GridSearchCV (sklearn) object for RandomForstClassification [1]

-score_testing_woparent_subsampling.csv
scoring result of [3]
-grid_cv_woparent_subsampling.pickle
GridSearchCV (sklearn) object for SVM [3]
-grid_KN_woparent_subsampling.pickle
GridSearchCV (sklearn) object for KNearestClassification [3]
-grid_RF_woparent_subsampling.pickle
GridSearchCV (sklearn) object for RandomForstClassification [3]

