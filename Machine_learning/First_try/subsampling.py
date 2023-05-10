from imblearn.under_sampling import RandomUnderSampler
from sklearn.neighbors import KNeighborsClassifier
import pickle
import os
import matplotlib
import matplotlib.pyplot as plt
import collections as cl
import csv
import pandas as pd
from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import sklearn.metrics as met
from sklearn.tree import export_graphviz
from IPython.display import Image
import graphviz
from sklearn.tree import export_graphviz
from IPython.display import Image
import graphviz

df=pd.read_csv("~/M1/get_nodes/node.csv")
interet=["A","G","C","U","PSU","5MC","OMG","5MU","2MG","OMC","A2M"]
scoring = ['f1_weighted', 'balanced_accuracy','precision_weighted','neg_log_loss']
df=df.loc[df['reponse'].isin(interet)]
print("REPARTITION PARENT",cl.Counter(df.loc[~df['reponse'].isin(['A','C','G','U'])]['parent']))
one_hot = pd.get_dummies(df['parent'])
df = df.drop('parent',axis = 1)
#df = df.join(one_hot)
mat=df.to_numpy()
X=mat[:,1:]
y=mat[:,0]
print(X[:5])
print(y)
rus = RandomUnderSampler(random_state=0)
X_resampled, y_resampled = rus.fit_resample(X, y)
print('Original dataset shape %s' % cl.Counter(y))
print('New dataset shape %s' % cl.Counter(y_resampled))

def evaluate(pred_list,test):
    res={}
    res['model']=["SVC","KNClassifier","RandomForst"]
    res['f1']=[]
    res['balanced_accuracy_score']=[]
    res['average_precision_score']=[]
    res['average_recall_score']=[] 
    #res['log_loss']=[]
    
    for pred in pred_list:
        res['f1'].append(met.f1_score(test,pred,average='weighted'))
        res['balanced_accuracy_score'].append(met.balanced_accuracy_score(test,pred))
        res['average_precision_score'].append(met.precision_score(test,pred,average='weighted'))
        res['average_recall_score'].append(met.recall_score(test,pred,average='weighted'))
	#res['log_loss'].append(met.log_loss(test,pred))
    data=pd.DataFrame(res)  
    data.to_csv("~/M1/score_testing_woparent_subsampling.csv",header=True,index=False)
    return None
def research_best(X,y):
    x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=42)

    #SVC
    svc = svm.SVC()
    parameters = {'probability':[True],'kernel':('poly','linear', 'rbf'), 'C':[1,2,4,6,8,10,12]}
    clf = GridSearchCV(svc, parameters,scoring=scoring,refit="f1_weighted")
    clf.fit(x_train,y_train)
    with open(os.path.join(".","grid_cv_woparent.pickle"), 'wb') as handle:
        pickle.dump(clf, handle)
    #PRED SVC
    svc_best=svm.SVC(**clf.best_params_)
    svc_best.fit(x_train,y_train)
    pred_svc=svc_best.predict(x_test)
    #EVALUATION
    
    print("FIN CSV")
    neigh = KNeighborsClassifier(n_neighbors=len(interet))
    parameters = {'weights':('uniform','distance'),'algorithm':('ball_tree', 'kd_tree', 'brute'),'leaf_size':range(400)[10:400:50]}
    clf2 = GridSearchCV(neigh, parameters,scoring=scoring,refit="f1_weighted")
    clf2.fit(x_train,y_train)
    with open(os.path.join(".","grid_KN_woparent.pickle"), 'wb') as handle:
        pickle.dump(clf2, handle)
    
    knn_best=KNeighborsClassifier(**clf2.best_params_)
    knn_best.fit(x_train,y_train)
    pred_knn=knn_best.predict(x_test)
    
    print("FIN KN")

    randfor = RandomForestClassifier(random_state=0)
    parameters = {"criterion":('gini', 'entropy', 'log_loss'),'min_samples_leaf':range(50)[1:50:5]}
    clf3 = GridSearchCV(randfor, parameters,scoring=scoring,refit="f1_weighted")
    clf3.fit(x_train,y_train)
    with open(os.path.join(".","grid_RF_woparent.pickle"), 'wb') as handle:
        pickle.dump(clf3, handle)
    
    randfor_best=RandomForestClassifier(**clf3.best_params_)
    randfor_best.fit(x_train,y_train)
    pred_randfor=randfor_best.predict(x_test)
    for i in range(1):
    	tree = randfor_best.estimators_[i]
    	dot_data = export_graphviz(tree,
                               feature_names=df.columns[1:],  
                               filled=True,   
                               impurity=False, 
                               proportion=True,
                               out_file="random_forst_woparent"+str(i)+".dot")
    print("FIN RANDOMFORST")

    evaluate([pred_svc,pred_knn,pred_randfor],y_test)
    return

research_best(X_resampled,y_resampled)
