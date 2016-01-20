function low_rank_rnd_vector(d,q,n,rho,gap,slope,method,options)
%%
options.slope=1;
options.rho=.1;
options.gap=.1;

slope=options.slope;
rho=options.rho;
gap=options.gap;
%%
d=128*128;
q=50;
n=1500;
rho=0.1;
gap=0.5;
slope=1;
%%
Vb = orth(normrnd(0,1,d,q)); %compute orthonormal basis
sigm = sqrt(rho)+sqrt(gap)+slope*[q:-1:1];

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
plot(LATENT)
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
