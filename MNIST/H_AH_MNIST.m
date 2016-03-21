clear all

fea=1; % if 1 uses orl otherwise MNIST
if fea    
    load('ORL_32x32');
    x1=fea';
    mu=mean(x1,2);
    x=bsxfun(@minus,x1,mu);
    
else
    x1 = loadMNISTImages('train-images-idx3-ubyte');
    mu=mean(x1,2);
    x=bsxfun(@minus,x1,mu);
    x=x(:,1:20000);
end
% eps_std=.1;
% x=bsxfun(@times,x,1./(eps_std+std(x,[],2)));
if 0
    sigma = cov(x');
    k=128;
    [U,S,V] = svds(sigma,k);
    epsilon=0;
    x = U * diag(1./sqrt(diag(S) + epsilon)) * U' * x;
end
fea=0;
%% batch PCA
q=64;
[d,T]=size(x);
d1=sqrt(d);
tic
[W1,Y1_tot]=pca(x','NumComponents',q);
toc
W1=W1';
P_q=W1'*W1;
%% show the effects of decreasing projection error
n_iter=1;
if fea
    %in this case initialize with PCA otherwise poor convergence
    n_init_pca=q;
    [M,W,Ysq]=run_H_AH_PCA(x(:,1:n_init_pca+1),q,n_iter,n_init_pca);
else
    [M,W,Ysq]=run_H_AH_PCA(x(:,1:3),q,n_iter,0);
end

Y_tot=(eye(q)+M)\W*x;
P_u=orth(W')*orth(W')';
subplot(4,4,1)
proj_err=norm(P_q-P_u,'fro')/norm(P_q,'fro');
imagesc(reshape(mu+P_u*x(:,end),[d1 d1]))
counter=1;
index=3;
title(['Samples:' num2str(index) ', E=' num2str(proj_err)])
axis off
axis image

for kk=2:15
    counter=counter+1;
    subplot(4,4,counter)
    if(fea)
        % in this case go through the dataset several times
        [M,W,Ysq]=run_H_AH_PCA(x(:,max(1,mod(index+1:index+kk^3,400))),q,n_iter,0,W,M,Ysq);
    else
        [M,W,Ysq]=run_H_AH_PCA(x(:,index+1:index+kk^3),q,n_iter,0,W,M,Ysq);
    end
    P_u=orth(W')*orth(W')';
    proj_err=norm(P_q-P_u,'fro')/norm(P_q,'fro');
    imagesc(reshape(mu+P_u*x(:,end),[d1 d1]))
    index=index+kk^3;
    title(['Samples:' num2str(index) ', E=' num2str(proj_err)])
    axis off
    axis image
    drawnow
end
subplot(4,4,counter+1)
imagesc(reshape(mu+P_q*x(:,end),[d1 d1]))
axis off
axis image
%% run simulation comparing algorithms


%%
T1=kmeans(Y1_tot,10,'Replicates',10);
%%
Y_tot=(eye(q)+M)\W*x;
T=kmeans(Y_tot',10,'Replicates',10);
%%

Y_tot=(eye(q)+M)\W*x;
[val,idx]=max(Y_tot,[],1);
for jj=1:10
    subplot(4,3,jj)
    imagesc(reshape(mean(x1(:,T==jj),2),[d1,d1]))
end
%%
lab=loadMNISTLabels('train-labels-idx1-ubyte');
MIhat = MutualInfo(T1-1,lab)
MIhat = MutualInfo(T-1,lab)
%%

for jj=0:9
    subplot(4,3,jj+1)
    imagesc(reshape(median(x(:,lab==jj),2),[d1,d1]))
end