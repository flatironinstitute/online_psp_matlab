# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:25:31 2016

@author: agiovann
"""

#%%
import numpy as np
import scipy
#%%
def run_OSM_PCA(X1,q, n_iter=50,init_weight_mult=np.pi):
    X=X1.copy()-np.mean(X1,0)
    M=np.random.random((q,q))*0
    T,d=np.shape(X)
    W=np.random.random((d,q))/np.repeat(np.sqrt(d),d)[:,np.newaxis]
    W=W.T
    Ysq=(init_weight_mult*np.linalg.norm(X[0,:],ord=2)**2/d)*(np.ones([np.shape(W)[0],1]))
#    
    Ysq=np.ones(q)/q
    Ysq=np.atleast_2d(Ysq).T
    
    for n in range(n_iter):
        print 'iteration:' + str(n)
        Ys=[]
        for i,x in enumerate(X):
            
            M,W,Y,Ysq = OSM_PCA (M,W,Ysq,x,q)
            Ys.append(Y)
        
    return M,W,np.squeeze(np.array(Ys))
#%%
def  OSM_PCA (M,W,Ysq,x,q, tol=1e-5, mat_iter=100, lambda_=0):

    """
    function implementing Hebbian Anti-hebbian networks for subspace learning
    
    
    Parameters
    -----------
    W: forward connection matrix
    M: lateral connection matrix
    Y: output space
    
    % options.update_method: method to perform the update of Y variable. 'ls' (least-square),
    % 'coord_desc' (coordinate descent), or 'mat_mult updated' (update all coordinated simultaneously)
    % options.tol: tolerance on convergence.
    % options.lambda: factor influencing decorrelation (0 no decorrelatoin applied, see NIPS 2015)
    % reference: Pehlevan et al, Neural Computation, 2015. Pehlevan et al, NIPS, 2015   
    % gamma=1./Ysq;
    """
    
    
    
    x=np.atleast_2d(x).T
    d=np.shape(W)[-1]

    # least square
    Y=scipy.linalg.solve((np.eye(q)+M),W.dot(x))        
    Ysq = Ysq + Y**2
    
    # Update weights
    Y_tmp=Y/Ysq
    Y_tmp_sq=Y**2/Ysq;

    # Update weights    
    W = W +   np.outer(Y_tmp,x.T) - W*Y_tmp_sq
    if np.isnan(np.sum(W)):
        W[np.isnan(W)] = 0;    

    M = M + (1+lambda_)*np.outer(Y_tmp,Y.T) - M*Y_tmp_sq   
    M[np.isnan(M)] = 0;
    
    # diagonal to zero
    M[::q+1]=0
    
    return M,W,Y,Ysq
    
#%%

#%%
if __name__ == "__main__":
    #%%
    n_components = 2
    import numpy as np
    import matplotlib.pyplot as plt
    
    from sklearn.datasets import load_iris
    from sklearn.decomposition import PCA, IncrementalPCA
    
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
            plt.title(title + " of iris dataset\nMean absolute unsigned error "
                      "%.6f" % err)
        elif "OSM" in title:
            err = np.abs(np.abs(X_pca) - np.abs(X_osmpca)).mean()
            plt.title(title + " of iris dataset\nMean absolute unsigned error "
                      "%.6f" % err)
        else:
            plt.title(title + " of iris dataset")
        plt.legend(loc="best")
        plt.axis([-4, 4, -1.5, 1.5])
    
        plt.show()
        
