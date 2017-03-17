# NeuralPCA
Benchmark of online PCA algorithms 

Set of functions efficiently implemented in MATLAB

## Installation

Clone the repository or unzip the source and add recursively folders from the root to the MATLAB path

Example

``` Matlab

%%Use ONLINE PCA with Hebbian Anti-Hebbian Learning rule (CEnzig et. al,2015 Neural Computation)
clear all
load fisheriris
x = bsxfun(@minus,meas,mean(meas,2));

n_init_PCA = 0; % number of samples used for initialization
d=size(x,2); % input dimensionality
q=2; % output dimensionality
n=size(x,1);
scrambled = randperm(n);
x = x(scrambled,:);
%% BATCH PCA

[coeff,score,pcvar] = pca(x,'NumComponents',2);

%% ONLINE PCA
[M,W,Ysq,Y]=run_H_AH_PCA(x',q,n_init_PCA,[],[],[]);
[M2,W2,Ysq2,Y2]=run_H_AH_PCA(x',q,n_init_PCA,W,M,Ysq);

% if you want to reconver the orthogonal projection
orth_projection = (pinv(diag(ones(q,1))+M(1:q,1:q))*W(1:q,:))';
[~,idx] = sort(sqrt(mean(1./Y.^2,2)),'descend');
%% PLOT RESULT
figure
scatterhist(score(:,1),score(:,2), 'group',species(scrambled))
%%
figure
scatterhist(Y(idx(1),:),Y(idx(2),:), 'group',species(scrambled(n_init_PCA+1:end)))
%%
figure
scatterhist(Y2(idx(1),:),Y2(idx(2),:), 'group',species(scrambled(n_init_PCA+1:end)))


```

## References

[1] Pehlevan, Cengiz, Tao Hu, and Dmitri B. Chklovskii. "A Hebbian/Anti-Hebbian Neural Network for Linear Subspace Learning: A Derivation from Multidimensional Scaling of Streaming Data." Neural computation (2015)

[2] Cardot, Hervé, and David Degras. "Online Principal Component Analysis in High Dimension: Which Algorithm to Choose?." arXiv preprint arXiv:1511.03688 (2015).
