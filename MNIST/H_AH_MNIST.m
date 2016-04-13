clear all
% close all
q=64;
fea=0; % if 1 uses orl otherwise MNIST
if fea    
    load('ORL_32x32');
%     load('YaleB_32x32');
%     load('ATT_faces_112_92');
%     d1=112;
%     d2=92;
%     frame rate is 120 Hz but frames might not be equally spaced
%     fea = hdf5read('agchr2_030915_01_040215_a_location_2_ds.hdf5','mov');
%     fea = ipermute(fea,[2 1 3]); % necessary since matlab and Python handle differntly the DHF5
 %    [d1,d2,T]=size(fea);
%     fea=reshape(fea,[d1*d2,T])';

    x1=fea';
%    x1=x1/max(abs(x1(:)));
    %
%     x=x1;
    mu=mean(x1,2);    
    x=bsxfun(@minus,x1,mu);
    
    
    %x=x/5/std(x(:));
%     x=x/range(x(:));
    
%     mu=mean(x,2);
%     x=x/quantile(abs(x(:)),.95);
    % OK!!
%      x=2*(x-min(x(:)))/(max(x(:))-min(x(:)))-1;   
%      x=x*sqrt(q);
%      mu=mean(x,2);
%      x=bsxfun(@minus,x,mu);
else
    x1 = loadMNISTImages('train-images-idx3-ubyte');
    x1=x1(:,1:20000);
   % x=x1;
   % mu=0;
%     x1=x1/max(abs(x1(:)));
     mu=mean(x1,2);        
     x=bsxfun(@minus,x1,mu);
     %x=x/3/std(x(:));
     
     %mu=mean(x,2);  
      
     x=2*(x-min(x(:)))/(max(x(:))-min(x(:)))-1;
     x=x*sqrt(q);
     mu=mean(x,2);
     x=bsxfun(@minus,x,mu);
end
% eps_std=.1;
% x=bsxfun(@times,x,1./(eps_std+std(x,[],2)));
if 0
    sigma = cov(x');
    k=q;
    [U,S,V] = svd(sigma);
    epsilon=0;
%     x = U * diag(1./sqrt(diag(S) + epsilon)) * U' * x;
    SS=diag(S);
    x = U * diag(1./[sqrt(SS(1:k) + epsilon)/SS(1); ones(960,1)]) * U' * x;
    x=x/max(abs(x(:)));
    mu=mean(x,2);
    x=bsxfun(@minus,x,mu);
end
total_img=size(x,2);
x=x(:,randperm(total_img));
img_test=x(:,end);
x=x(:,1:end-1);
total_img=total_img-1;
%% batch PCA
[d,T]=size(x);
if ~exist('d1')
    d1=sqrt(d);
    d2=sqrt(d);
end

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
%%
tic
rec_err_real=norm(x-P_q*x,'fro')^2/norm(x,'fro')^2
toc
%%
tic
nsamp=50;
errs=zeros(1,nsamp);
idx=randsample(size(x,2),nsamp);
for kk=1:nsamp
    errs(kk)=norm(x(:,idx(kk))-P_q*x(:,idx(kk)))^2/norm(x(:,idx(kk)))^2;
%     errs(kk)=norm(x(:,kk)-P_q*x(:,kk))^2/norm(mu+x(:,kk))^2;
end
rec_err_real=mean(errs)
toc
%% show the effects of decreasing projection error
figure
n_iter=1;
init_iter=q;
is_method_OSM=1; %OSM or IPCA
loop_data=1;
init_pca=0;
tic

if init_pca
    %in this case initialize with PCA otherwise poor convergence
    if is_method_OSM==1
        [M,W,Ysq]=run_H_AH_PCA(x(:,1:init_iter+1),q,init_iter,[],[],[]);
        
    else
        [values,vectors,iter_so_far]=run_incrPCA(x(:,1:init_iter+1),q,init_iter,[],[],[]);
    end
else
    if is_method_OSM==1
        [M,W,Ysq]=run_H_AH_PCA(x(:,1:init_iter),q,0,[],[],[]);
    else
        [values,vectors,iter_so_far]=run_incrPCA(x(:,1:init_iter),q,0,[],[],[]);
    end
end

% Y_tot=(eye(q)+M)\W*x;
if is_method_OSM == 1
      F=(pinv(eye(q)+M)*W); 
      
      P_u=orth(F')*orth(F')';  
      proj_err=compute_projection_error(W1',orth(F'));
else
      proj_err=compute_projection_error(W1',vectors);
      P_u=vectors*vectors';
end

subplot(4,4,1)
% proj_err=norm(P_q-P_u,'fro')/norm(P_q,'fro');

imagesc(reshape(mu+img_test,[d1 d2]))
%hist(M(:),100)
counter=1;
index=init_iter;
title(['S:' num2str(index) ', E:' num2str(proj_err)])
 axis off
 axis image

for kk=2:15
    counter=counter+1;    
    subplot(4,4,counter)
    if(loop_data)
        % in this case go through the dataset several times
        if is_method_OSM==1
             [M,W,Ysq]=run_H_AH_PCA(x(:,max(1,mod(index+1:index+kk^3,total_img))),q,0,W,M,Ysq);
        else
             [values,vectors,iter_final]=run_incrPCA(x(:,max(1,mod(index+1:index+kk^3,total_img))),q,0,values,vectors,index);
        end
        upper=index+kk^3;
    else
        if index>total_img
            break
        end
        upper=min(total_img,index+kk^3);
        if is_method_OSM == 1 
            [M,W,Ysq]=run_H_AH_PCA(x(:,index+1:upper),q,0,W,M,Ysq);
        else
            [values,vectors,iter_final]=run_incrPCA(x(:,index+1:upper),q,0,values,vectors,index);
        end
    end
    
   
   if is_method_OSM == 1
        F=(pinv(eye(q)+M)*W);
        P_u=orth(F')*orth(F')'; 
        rec_err=norm(x-P_u*x,'fro')^2/norm(x,'fro')^2          
        proj_err=compute_projection_error(W1',orth(F'))
   else
        proj_err=compute_projection_error(W1',vectors)        
        P_u=vectors*vectors';
        rec_err=norm(x-P_u*x,'fro')^2/norm(x,'fro')^2;
   end
   
%     errs=zeros(1,nsamp);
%     idx=randsample(size(x,2),nsamp);
%     for kkk=1:nsamp
%         errs(kkk)=norm(x(:,idx(kkk))-P_q*x(:,idx(kkk)))^2/norm(mu+x(:,idx(kkk)))^2;
%     end
%     rec_err=mean(errs);

%     proj_err=norm(P_q-P_u,'fro')/norm(P_q,'fro');
    imagesc(reshape(mu+P_u*img_test,[d1 d2]))
%    hist(M(:),100) 
    index=index+kk^3;
    % round(1000*rec_err/rec_err_real)/1000
    title(['S:' num2str(round(upper/10)) ', E:' num2str(round(1000*proj_err)/1000) ', R:' num2str(rec_err)])
     axis off
     axis image
     colormap gray
    drawnow
end

subplot(4,4,counter+1)
imagesc(reshape(mu+P_q*img_test,[d1 d2]))
axis off
axis image
toc
%% run simulation comparing algorithmsnsamp
nsamp=100;
errs=zeros(1,nsamp);
for kk=1:nsamp
    errs(kk)=norm(x(:,kk)-P_u*x(:,kk))^2/norm(mu+x(:,kk))^2;
end
rec_err_real=mean(errs)

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