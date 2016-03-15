clear all
x1 = loadMNISTImages('train-images-idx3-ubyte');
x=bsxfun(@minus,x1,mean(x1,2));
% eps_std=.1;
% x=bsxfun(@times,x,1./(eps_std+std(x,[],2)));
if 0
    sigma = x * x' / size(x, 2);
    k=128;
    [U,S,V] = svds(sigma,k);
    epsilon=0;
    x = U * diag(1./sqrt(diag(S) + epsilon)) * U' * x;
end
%rec_error=norm(x-U*U'*x,2)/norm(x,2)

%%
tic
[d,T]=size(x);
q=64;
W=rand(q,d)*.01;
M=rand(q,q)*.01;
M=M-diag(diag(M));

Y=rand(q,1)*0;

options_algorithm=struct();
options_algorithm.pca_algorithm='H_AH_NN_PCA';
options_algorithm.q=q;
options_algorithm.update_method='ls';
options_algorithm.tol=1e-5;

Ysq=10*ones(size(W,1),1);

for iter=1:1
    disp(iter)
    scramble=randperm(T);
    for kk=1:T
        options_algorithm.gamma=1./Ysq;
        [M,W,Y]=H_AH_NN_PCAFast(M,W,Y,x(:,scramble(kk))',options_algorithm);
        Ysq = Ysq + Y.^2;
    end
end
toc
%%
tic
[W1,Y1_tot]=pca(x','NumComponents',q);
% [Y1_tot,W1]=nnmf(x',q);

W1=W1';
toc
%%
norm(x-W1'*W1*x)/norm(x)
%%
norm(x-orth(W')*orth(W')'*x)/norm(x)
%%
T1=kmeans(Y1_tot,10,'Replicates',10);
%%
Y_tot=(eye(q)+M)\W*x;
T=kmeans(Y_tot',10,'Replicates',10);
%%
d1=sqrt(d);
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