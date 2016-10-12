%%Use ONLINE PCA with Hebbian Anti-Hebbian Learning rule (CEnzig et. al,2015 Neural Computation)
clear all

tic
n0 = 10; % number of samples used for initialization
d=1000; % input dimensionality
q=200; % output dimensionality
n=1000;

% create covariance matrix for brownian motion

% hold all;
[I,J] = ind2sub([d,d],1:d^2);
C=zeros(d);
C(:)=min(I,J)/d; % covariance matrix
[Q,L,~]=svd(C); 
Pq=Q(:,1:q)*Q(:,1:q)'; % projection vector

x1= normrnd(0,9/sqrt(d),n,d);
x= cumsum(x1,2);

% initialize with batch PCA
disp('svd')
[V,D,W] = svd(cov(x(1:n0,:)));
%eigval=diag(D);
values=diag(D);
eigvect=V;
%%
profile on
tic
values=values(1:q);
vectors=eigvect(:,1:q);



M = 0*randn(q,q); 
W = vectors;
W=W';
Y = zeros(q,1);
Ysq=sum((W*x(1:n0,:)').^2,2);% learning rate

for ii = 1:q
    M(ii,ii) = 0;
end

errors=zeros(n,1)*NaN;
tol=1e-5;
tms=zeros(1,n-n0);
% tms=[];
options.tol=1e-5;

[eig_vect,~,eig_val]=pca(x,'NumComponents',q);

for i = (n0+1):n
    if mod(i,100)==1
        disp(i)
    end    
    tic
    options.q=q;
    options.n=i-1;    
    [values, vectors] = CCIPCA(values, vectors, x(i,:),options);
%     [values, vectors] = incrPCA_fast(values, vectors, x(i,:),options);    
    tms(i-n0)=toc;
    errors(i,:)=compute_projection_error(eig_vect,orth(vectors));

end
hold all
plot(errors)

ylabel('Eigenspace Estimation Error Pq vs Pqhat')
xlabel('Iterations')
toc
profile off
profile viewer
% figure
% hist(tms,0:.0001:.01),xlim([0 .01])
% ylabel('Counts')
% xlabel('Time (s)')