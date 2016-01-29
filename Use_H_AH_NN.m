%%Use ONLINE PCA with Hebbian Anti-Hebbian Learning rule (CEnzig et. al,2015 Neural Computation)
clear all

tic
n0 = 50; % number of samples used for initialization
d=5000; % input dimensionality
q=50; % output dimensionality
n=1000;

% create covariance matrix for brownian motion

hold all;
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

errors=zeros(n,1);
tol=1e-5;
tms=zeros(1,n-n0);
% tms=[];
options.tol=1e-5;
options.update_method='ls';


for i = (n0+1):n
%     disp(i)
    tic
    %[values, vectors] = incrPCA_fast(values, vectors, x(i,:),i-1, [], q);    
     [M,W,Y]=H_AH_NN_fast(x(i,:)',Y,W,M,q,Ysq,options);
     Ysq = Ysq + Y.^2; % update learning rate
    tms(i-n0)=toc;
%      if mod(i,50)==0
%          tms=[tms toc];
%          tic
%      end
%     if mod(i,100)==1
%        vectors = (pinv(diag(ones(q,1))+M(1:q,1:q))*W(1:q,:))';
%        Pq_hat=vectors(:,1:q)*vectors(:,1:q)';
%        errors(i)=(norm(Pq_hat-Pq,'fro')^2)/(norm(Pq,'fro')^2);
%     end
    %         errors(i,ll)=2*(1-trace(Pq_hat*Pq)/q);
end
hold all
plot(errors)

ylabel('Eigenspace Estimation Error Pq vs Pqhat')
xlabel('Iterations')
toc
profile off
profile viewer
figure
hist(tms,0:.0001:.01),xlim([0 .01])
ylabel('Counts')
xlabel('Time (s)')