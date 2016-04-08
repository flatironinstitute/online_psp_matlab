clear all

fea=0; % if 1 uses orl otherwise MNIST
if fea    
    load('ORL_32x32');
    x1=fea';
    x1=x1/max(abs(x1(:)));
    mu=mean(x1,2);
    x=bsxfun(@minus,x1,mu);
%     x=2*(x-min(x(:)))/(max(x(:))-min(x(:)))-1;    
else
    x1 = loadMNISTImages('train-images-idx3-ubyte');
    x1=x1(:,1:20000);
   % x=x1;
   % mu=0;
%     x1=x1/max(abs(x1(:)));
     mu=mean(x1,2);        
     x=bsxfun(@minus,x1,mu);
        
end
% eps_std=.1;
% x=bsxfun(@times,x,1./(eps_std+std(x,[],2)));
if 0
    sigma = cov(x');
    k=64;
    [U,S,V] = svd(sigma);
    epsilon=0;
%     x = U * diag(1./sqrt(diag(S) + epsilon)) * U' * x;
    SS=diag(S);
    x = U * diag(1./[sqrt(SS(1:k) + epsilon)/SS(1); ones(960,1)]) * U' * x;

end

%% batch PCA
q=64;
[d,T]=size(x);
d1=sqrt(d);

% evl=eig(cov(x'));
% evl(end:end-5)
% 
% [S,V,D]=svd(cov(x'));
% sevl=diag(V);
% sevl(1:5)

tic
[W1,Y1_tot,evl]=pca(x','NumComponents',q);
evl=evl(1:q);
toc
W1=W1';
P_q=W1'*W1;

%% show the effects of decreasing projection error
figure
n_iter=1;
init_iter=q;
% method=0; %IPCA
method=0; %OSM
fea=0;
if fea
    %in this case initialize with PCA otherwise poor convergence

    if method==1
        [M,W,Ysq]=run_H_AH_PCA(x(:,1:init_iter+1),q,init_iter);
    else
        [values,vectors,iter_so_far]=run_incrPCA(x(:,1:init_iter+1),q,init_iter,[],[],[]);
    end
else
    if method==1
        [M,W,Ysq]=run_H_AH_PCA(x(:,1:init_iter),q,0);
    else
        [values,vectors,iter_so_far]=run_incrPCA(x(:,1:init_iter),q,0,[],[],[]);
    end
end

% Y_tot=(eye(q)+M)\W*x;
if method == 1
      F=(pinv(eye(q)+M)*W);
      P_u=orth(F')*orth(F')';  
     proj_err=compute_reconstruction_error(W1',orth(F'))
else
    proj_err=compute_reconstruction_error(W1',vectors)
     P_u=vectors*vectors';
end

subplot(6,5,1)
% proj_err=norm(P_q-P_u,'fro')/norm(P_q,'fro');

imagesc(reshape(mu+P_u*x(:,end),[d1 d1]))
counter=1;
index=init_iter;
title(['Samples:' num2str(index) ', E=' num2str(proj_err)])
axis off
axis image

for kk=2:16
    counter=counter+1;
    
    subplot(6,5,counter)
    fea=1;
    if(fea)
        % in this case go through the dataset several times
        if method==1
             [M,W,Ysq]=run_H_AH_PCA(x(:,max(1,mod(index+1:index+kk^3,400))),q,0,W,M,Ysq);
        else
             [values,vectors,iter_final]=run_incrPCA(x(:,max(1,mod(index+1:index+kk^3,400))),q,0,values,vectors,index);
        end
    else
        if method == 1 
            [M,W,Ysq]=run_H_AH_PCA(x(:,index+1:index+kk^3),q,0,W,M,Ysq);
        else
            [values,vectors,iter_final]=run_incrPCA(x(:,index+1:index+kk^3),q,0,values,vectors,index);
        end
    end
    
   if method == 1
        F=(pinv(eye(q)+M)*W);
      P_u=orth(F')*orth(F')';  
        proj_err=compute_reconstruction_error(W1',orth(F'))
   else
        proj_err=compute_reconstruction_error(W1',vectors)
     P_u=vectors*vectors';
    end

%     proj_err=norm(P_q-P_u,'fro')/norm(P_q,'fro');
    imagesc(reshape(mu+P_u*x(:,end),[d1 d1]))
%    hist(M(:),100)
    index=index+kk^3;
    title(['Samples:' num2str(index) ', E=' num2str(proj_err)])
     axis off
     axis image
    drawnow
end
subplot(6,5,counter+1)
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