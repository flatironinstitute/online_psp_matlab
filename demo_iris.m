%%Use ONLINE PCA with Hebbian Anti-Hebbian Learning rule (CEnzig et. al,2015 Neural Computation)
clear all
load fisheriris
x = bsxfun(@minus,meas,mean(meas,1));

n_init_PCA = 0; % number of samples used for initialization
d=size(x,2); % input dimensionality
q=2; % output dimensionality
n=size(x,1);
scrambled = randperm(n);
x = x(scrambled,:);
%% BATCH PCA

[coeff,score,pcvar] = pca(x,'NumComponents',2);

%% ONLINE PCA
[M,W,Ysq,Y]=run_OSM_PCA(x',q,n_init_PCA,[],[],[]);
% pass over data once more to evaluate improvement
[M2,W2,Ysq2,Y2]=run_OSM_PCA(x',q,0,W,M,Ysq);

% if you want to reconver the orthogonal projection
% orth_projection = (pinv(diag(ones(q,1))+M(1:q,1:q))*W(1:q,:))';
%% PLOT RESULT BATCH
figure
scatterhist(score(:,1),score(:,2), 'group',species(scrambled))
%% FIRST ITERATION
figure
scatterhist(Y(1,:),Y(2,:), 'group',species(scrambled(n_init_PCA+1:end)))
%% SECOND ITERATION
figure
scatterhist(Y2(1,:),Y2(2,:), 'group',species(scrambled(n_init_PCA+1:end)))