#coding=utf-8
#!/usr/bin/env python3

import sys
import os
import random
#import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import ROOT
import joblib

import imblearn
from imblearn.combine import SMOTETomek
from imblearn.over_sampling import SMOTE

from sklearn import datasets
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import classification_report, roc_auc_score,accuracy_score
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split


import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import learning_curve
from sklearn.model_selection import ShuffleSplit
from sklearn.svm import SVC


import xgboost as xgb

def usage():
	print ('test usage')
	sys.stdout.write('''
			SYNOPSIS

			./BDT_pre.py signal bkg 

			AUTHOR
			Yanxi Gu <GUYANXI@ustc.edu>

			DATE
			10 Jan 2021
			\n''')

def main():

    args = sys.argv[1:]
    if len(args) < 2:
        return usage()

    print ('part1')   

    # get root files and convert them to array
    #branch_names = """Px_Z,Py_Z,Pz_Z,E_Z,Px_H,Py_H,Pz_H,E_H,Px_H_Z,Py_H_Z,Pz_H_Z,E_H_Z,Px_H_Zs,Py_H_Zs,Pz_H_Zs,E_H_Zs,Px_Z_Mup,Py_Z_Mup,Pz_Z_Mup,E_Z_Mup,Px_Z_Mum,Py_Z_Mum,Pz_Z_Mum,E_Z_Mum,Px_H_Z_Mup,Py_H_Z_Mup,Pz_H_Z_Mup,E_H_Z_Mup,Px_H_Z_Mum,Py_H_Z_Mum,Pz_H_Z_Mum,E_H_Z_Mum""".split(",")
    #branch_names = """Px_Z,Py_Z,Pz_Z,E_Z,Px_H,Py_H,Pz_H,E_H""".split(",")
    #branch_names = """costheta1,costheta2,phi,M_H,M_Z,M_H_Z,M_H_Zs,M_Z_Mup,M_Z_Mum""".split(",")
    #branch_names = """costheta1,costheta2,phi,phi1,costheta1_H,costheta2_H,phi_H""".split(",")
    #branch_names = """costheta1,costheta2,phi,P_H_Z,P_H_Zs,P_Z_Mup,P_Z_Mum,Px_Z,Py_Z,Pz_Z,Px_H,Py_H,Pz_H""".split(",")  #The last selected feature for training the truth value
    branch_names = """costheta1,costheta2,Px_H,Py_H,Pz_H,Px_Z,Py_Z,Pz_Z,Px_Z_Mup,Py_Z_Mup,Pz_Z_Mup,Pz_Z_Mup,E_Z_Mup,Px_Z_Mum,Py_Z_Mum,Pz_Z_Mum,E_Z_Mum""".split(",") # new sample
#    branch_names = """Px_Beamp,Py_Beamp,Pz_Beamp,E_Beamp,Px_Beamm,Py_Beamm,Pz_Beamm,E_Beamm,Px_Z,Py_Z,Pz_Z,E_Z,Px_H,Py_H,Pz_H,E_H,Px_H_Z,Py_H_Z,Pz_H_Z,E_H_Z,Px_H_Zs,Py_H_Zs,Pz_H_Zs,E_H_Zs,Px_Z_Mup,Py_Z_Mup,Pz_Z_Mup,E_Z_Mup,Px_Z_Mum,Py_Z_Mum,Pz_Z_Mum,E_Z_Mum,Px_H_Z_Mup,Py_H_Z_Mup,Pz_H_Z_Mup,E_H_Z_Mup,Px_H_Z_Mum,Py_H_Z_Mum,Pz_H_Z_Mum,E_H_Z_Mum""".split(",")

    fin1 = ROOT.TFile(args[0])
    fin2 = ROOT.TFile(args[1])

    tree1 = fin1.Get("trialTree")   #truth's root tree
    #tree1 = fin1.Get("fancy_tree") #Reconstruction's root  tree
    signal0 = tree1.AsMatrix(columns=branch_names)
    signal = signal0[:100000,:]
    #signal = signal0[:100000,:]
    tree2 = fin2.Get("trialTree") #truth's root tree
    #tree2 = fin2.Get("fancy_tree") #Reconstruction's root  tree
    backgr0 = tree2.AsMatrix(columns=branch_names)
    backgr = backgr0[:100000,:]
    #backgr = backgr0[:100000,:]

    print('signal')
    print(signal)
    print('backgr')
    print(backgr)

    # for sklearn data is usually organised into one 2D array of shape (n_samples * n_features)
    # containing all the data and one array of categories of length n_samples
    X_raw = np.concatenate((signal, backgr))
    y_raw = np.concatenate((np.ones(signal.shape[0]), np.zeros(backgr.shape[0])))
    print(len(signal))
    print(len(backgr))

    print ('part2')

    #imbalanced learn
    n_sig = len(y_raw[y_raw==1])
    n_bkg = len(y_raw[y_raw==0])
    print(n_sig)
    print(n_bkg)
    sb_ratio = len(y_raw[y_raw==1])/(1.0*len(y_raw[y_raw==0]))
    if (sb_ratio > 0.2 and sb_ratio < 0.5):
        smote = SMOTE(ratio=0.5)
        X, y = smote.fit_sample(X_raw, y_raw)
        print ('Number of events: ')
        print ('before: signal: ', len(y_raw[y_raw==1]), ' background: ', len(y_raw[y_raw==0]))
        print ('after: signal: ', len(y[y==1]), ' background: ', len(y[y==0]))
    elif (n_sig > 1000 and sb_ratio < 0.2 and sb_ratio > 0.1):
        smote = SMOTE(ratio=0.2)
        X, y = smote.fit_sample(X_raw, y_raw)
        print ('Number of events: ')
        print ('before: signal: ', len(y_raw[y_raw==1]), ' background: ', len(y_raw[y_raw==0]))
        print ('after: signal: ', len(y[y==1]), ' background: ', len(y[y==0]))
    elif (n_sig < 1000 and sb_ratio < 0.2 and sb_ratio > 0.1):
        smote = SMOTE(ratio=0.4)
        X, y = smote.fit_sample(X_raw, y_raw)
        print ('Number of events: ')
        print ('before: signal: ', len(y_raw[y_raw==1]), ' background: ', len(y_raw[y_raw==0]))
        print ('after: signal: ', len(y[y==1]), ' background: ', len(y[y==0]))
    elif (sb_ratio < 0.1 and sb_ratio > 0.05):
        smote = SMOTE(ratio=0.4)
        X, y = smote.fit_sample(X_raw, y_raw)
        print ('Number of events: ')
        print ('before: signal: ', len(y_raw[y_raw==1]), ' background: ', len(y_raw[y_raw==0]))
        print ('after: signal: ', len(y[y==1]), ' background: ', len(y[y==0]))
    elif (sb_ratio < 0.05 and sb_ratio > 0.01):
        smote = SMOTE(ratio=0.1)
        X, y = smote.fit_sample(X_raw, y_raw)
        print ('Number of events: ')
        print ('before: signal: ', len(y_raw[y_raw==1]), ' background: ', len(y_raw[y_raw==0]))
        print ('after: signal: ', len(y[y==1]), ' background: ', len(y[y==0]))
    elif (sb_ratio < 0.01):
        smote = SMOTE(ratio=0.03)
        X, y = smote.fit_sample(X_raw, y_raw)
        print ('Number of events: ')
        print ('before: signal: ', len(y_raw[y_raw==1]), ' background: ', len(y_raw[y_raw==0]))
        print ('after: signal: ', len(y[y==1]), ' background: ', len(y[y==0]))
    else:
        X = X_raw
        y = y_raw
        print ('Number of events: ')
        print ('signal: ', len(y[y==1]), ' background: ', len(y[y==0]))


    """
    Training Part
    """
    # Train and test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.40, random_state=42)

    #dt = DecisionTreeClassifier(max_depth=51, min_samples_leaf=20, min_samples_split=40)

    #bdt = AdaBoostClassifier(dt, algorithm='SAMME', n_estimators=250, learning_rate=0.03)
    dt = DecisionTreeClassifier(max_depth=5, min_samples_leaf=100, min_samples_split=10)
   
    bdt = AdaBoostClassifier(dt, algorithm='SAMME', n_estimators=200, learning_rate=0.2)
    bdt.fit(X_train, y_train)

    importances = bdt.feature_importances_
    f = open('bdt_results/output_importance_New.txt', 'w')
    f.write("%-25s%-15s\n"%('Variable Name','Output Importance'))
    #for i in range (32):
    for i in range (17):
        f.write("%-25s%-15s\n"%(branch_names[i], importances[i]))
        print("%-25s%-15s\n"%(branch_names[i], importances[i]), file=f)
    f.close() 

    y_predicted = bdt.predict(X_train)
    print (classification_report(y_train, y_predicted, target_names=["background", "signal"]))
    print ("Area under ROC curve: %.4f"%(roc_auc_score(y_train, bdt.decision_function(X_train))))
    y_trainacc = accuracy_score(y_train, y_predicted)
    print("Area under ACC curve: %.4f"%y_trainacc)

    y_predicted = bdt.predict(X_test)
    print (classification_report(y_test, y_predicted, target_names=["background", "signal"]))
    print ("Area under ROC curve: %.4f"%(roc_auc_score(y_test, bdt.decision_function(X_test))))
    y_trainacc = accuracy_score(y_test, y_predicted)
    print("Area under ACC curve: %.4f"%y_trainacc)

    decisions1 = bdt.decision_function(X_train)
    decisions2 = bdt.decision_function(X_test)

    filepath = 'SM-vs-BSM-CPeven'

    # Compute ROC curve and area under the curve
    fpr1, tpr1, thresholds1 = roc_curve(y_train, decisions1)
    fpr2, tpr2, thresholds2 = roc_curve(y_test, decisions2)
    roc_auc1 = auc(fpr1, tpr1)
    roc_auc2 = auc(fpr2, tpr2)
    fig=plt.figure(figsize=(8,6))
    fig.patch.set_color('white')
    plt.plot(fpr1, tpr1, lw=1.2, label='train:ROC (area = %0.4f)'%(roc_auc1), color="r")
    plt.plot(fpr2, tpr2, lw=1.2, label='test: ROC (area = %0.4f)'%(roc_auc2), color="b")
    plt.plot([0,1], [0,1], '--', color=(0.6, 0.6, 0.6), label = 'Luck')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating  characteristic')
    plt.legend(loc = "lower right")
    plt.grid()
    plt.savefig('./bdt_results/'+filepath+'/ROC_Hbb.png')
#    plt.show()

    compare_train_test(bdt, X_train, y_train, X_test, y_test, filepath)

    joblib.dump(bdt, './bdt_results/'+filepath+'/bdt_model_New.pkl')

# Comparing train and test results
def compare_train_test(clf, X_train, y_train, X_test, y_test, savepath, bins=30):

    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)):
        d1 = clf.decision_function(X[y>0.5]).ravel()
        d2 = clf.decision_function(X[y<0.5]).ravel()
        decisions += [d1, d2]

    low = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    low_high = (low, high)
    fig=plt.figure(figsize=(8,5.5))
    fig.patch.set_color('white')
    plt.hist(decisions[0], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density = True, label='Signal (train)')
    plt.hist(decisions[1], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density = True, label='Background (train)')
    
    hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True)
    scale = len(decisions[2])/sum(hist)
    err = np.sqrt(hist*scale)/scale

    width = (bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='Signal (test)')

    hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
    scale = len(decisions[2])/sum(hist)
    err = np.sqrt(hist*scale)/scale

    plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='Background (test)')
  
    plt.xlabel("BDT score")
    plt.ylabel("Normalized Unit")
    plt.legend(loc='best')
    plt.savefig("./bdt_results/"+savepath+"/BDTscore_Hbb.png")
#    plt.show()

if __name__ == '__main__':
	print('start')
	main()
