function low_rank_rnd_vector(d,q,rho,gap)
%%
clear all
d=512*512;
q=50;
n=1500;
rho=0.1;
gap=0.5;
Vb = orth(normrnd(0,1,d,q)); %compute orthonormal basis
sigm = ones(q,1)*sqrt(gap);
%sigm(1:q)=sqrt([1 .9 .8 .7 .6]);
X=zeros(d,n);
disp('created X')
x = normrnd(0,1,q,n);
disp('created x')
eta = normrnd(0,1,d,n);
disp('created eta')
for kk=1:n
    if mod(kk,100)==0
        disp(kk)
    end
    X(:,kk)=sum(Vb*diag(sigm.*x(:,kk)),2)+sqrt(rho)*eta(:,kk);
end
%%
% C=cov(X');
% disp('created C')
% 
% [vectors,lam] = eig(C,'vector');
% disp('EIG')
% 

% %vectors=vectors(end:-1:end-4,:)';
% vectors=vectors(:,end:-1:end-4);
% [U,S,V] = svd(C,'econ');
% disp('SVD')
[COEFF, SCORE, LATENT] = pca(X','NumComponents',q);
disp('PCA')

%%
close all
subplot(2,2,1)
hold all
plot(sigm.^2+rho)
plot(LATENT)
subplot(2,2,2)
imagesc(Vb*Vb')
% subplot(2,3,3)
% imagesc(vectors*vectors',[-.3 .3])
subplot(2,2,3)
imagesc(COEFF*COEFF')
subplot(2,2,4)
imagesc(abs(COEFF*COEFF'-Vb*Vb'),[0 .01])
