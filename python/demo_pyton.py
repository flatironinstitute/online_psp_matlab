#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 11:37:37 2017

@author: agiovann
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_iris
from sklearn.decomposition import PCA, IncrementalPCA
from OSM import run_OSM_PCA
#%%
n_components = 2

iris = load_iris()
X = iris.data
y = iris.target


ipca = IncrementalPCA(n_components=n_components, batch_size=10)
X_ipca = ipca.fit_transform(X)
#%%
pca = PCA(n_components=n_components)
X_pca = pca.fit_transform(X)
#%%
M,W,Ys  = run_OSM_PCA(X,n_components)

F=(np.linalg.pinv(np.eye(n_components)+M[:n_components ,:n_components ]).dot(W[:n_components ,:])).T
X_osmpca = (X-np.mean(X,0)).dot(F)
X_osmpca = Ys
#%%
plt.close('all')
for X_transformed, title in [(X_osmpca, "OSM PCA"), (X_ipca, "Incremental PCA"), (X_pca, "PCA")]:
    plt.figure(figsize=(8, 8))
    for c, i, target_name in zip("rgb", [0, 1, 2], iris.target_names):
        plt.scatter(X_transformed[y == i, 0], X_transformed[y == i, 1],
                    c=c, label=target_name)

    if "Incremental" in title:
        err = np.abs(np.abs(X_pca) - np.abs(X_ipca)).mean()
#        plt.title(title + " of iris dataset\nMean absolute unsigned error "
#                  "%.6f" % err)
    elif "OSM" in title:
        err = np.abs(np.abs(X_pca) - np.abs(X_osmpca)).mean()
#        plt.title(title + " of iris dataset\nMean absolute unsigned error "
#                  "%.6f" % err)
    
    plt.title(title + " of iris dataset")
    plt.legend(loc="best")
    plt.axis([-4, 4, -1.5, 1.5])

    plt.show()
    
