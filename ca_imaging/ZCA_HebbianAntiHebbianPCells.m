%%
clear all
load('mov.mat')
%%
clear all
nam = 'demoMovie.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=2000;					% user input: how many frames to read   (optional, default until the end)

mov = bigread2(nam,sframe,num2read);
if ~isa(mov,'double');    mov = double(mov);  end         % convert to double
mov=mov(20:35,20:35,:);
%%

[d1,d2,T]=size(mov);

M=reshape(mov,[d1*d2 T])';
% imagesc(Yr./max(Yr(:)))
M(isnan(M))=0;
%[wcoeff,score,latent,tsquared,explained]  = princomp(M);
M=bsxfun(@minus,M,mean(M,1));
x=M;
%%
tic
sigma = x * x' / size(x, 2);
k=30;
[U,S,V] = svds(sigma,k);

xRot = U' * x;          % rotated version of the data.

xTilde = U(:,1:k)' * x; % reduced dimension representation of the data,
% where k is the number of eigenvectors to keep
epsilon=n0;
zca_whitesig = U * diag(1./sqrt(diag(S) + epsilon)) * U' * x;
% zca_whitesig = diag(1./sqrt(diag(S) + epsilon)) * U' * x;


toc
rec_error=norm(x-U*U'*x,2)/norm(x,2)
%%
tic
k1=3
[U,S,V] = svds(M,k1);
% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
whiteningMatrix = inv (sqrt(S)) * U';
dewhiteningMatrix = U * sqrt(S);
% Project to the eigenvectors of the covariance matrix.
% Whiten the samples and reduce dimension simultaneously.
whitesig =  whiteningMatrix * M;
toc
rec_error=norm(M-U*U'*M,2)/norm(M,2)

%% calculate show masks
for kk=1:k1
    subplot(2,1,1)
    %     zca_wsigmask=reshape(zca_whitesig',[d1 d2 30]);
    %     imagesc(zca_wsigmask(:,:,kk))
    axis image
    subplot(2,1,2)
    wsigmask=reshape(whitesig',[d1 d2 k1]);
    imagesc(wsigmask(:,:,kk))
    axis image
    colormap gray
    pause
end
%%

orgw=whiteningMatrix'*whitesig;
org=dewhiteningMatrix*whitesig;

%%
[o1, o2,o3]=fastica(whitesig,'numOfIC',3);
% [zo1, zo2,zo3]=fastica(zca_whitesig,'numOfIC',30);
o1mask=reshape(o1',[d1 d2 3]);
% zo1mask=reshape(zo1',[d1 d2 30]);

for kk=1:30
    subplot(2,1,1)
    imagesc(o1mask(:,:,kk));
    axis image
    subplot(2,1,2)
    % imagesc(zo1mask(:,:,kk));
    axis image
    
    colormap gray
    pause
end
%%

q=5;
d=d1*d2;
W=rand(q,d)*.1;
M=rand(q,q)*.1;
Y=rand(q,1)*.1;
options_algorithm=struct();
options_algorithm.pca_algorithm='H_AH_NN_PCA';
options_algorithm.q=q;
options_algorithm.update_method='ls';
options_algorithm.tol=1e-5;
options_algorithm.lambda=.5;


Ysq=10*ones(size(W,1),1);


for iter=1:10
    disp(iter)
    scramble=randperm(T);
    
    for kk=1:T
        options_algorithm.gamma=1./Ysq;
        [M,W,Y]=H_AH_NN_PCAFast(M,W,Y,x(scramble(kk),:),options_algorithm);
        Ysq = Ysq + Y.^2;
    end
end

for kk=1:size(W,1)
    mask=reshape(W(kk,:),[d1 d2]);
    imagesc(mask)
    axis image
    pause
end
%%
close all
%X=x';
X=bsxfun(@times,x',1./(std(x',[],2)));
lambda=500000; % for patch
% X=(x')./(max(x(:))-min(x(:)));
% lambda=1.5;
%X=X/norm(X,'fro');
% X=zca_whitesig';
% lambda=250;
clc
tic
T=size(X,2);
m=10;
n=size(X,1);

w_nonneg=1;
tolerance=.0001;
active_set=[];
W=zeros(m,n);
M=zeros(m,m);
y_t=zeros(m,1);
y_t_hat=zeros(m,1);
detected_clusters=[];
newX=[];
y_ts=[];
pixels=randperm(T);
pixels=1:T;
%
for nnn=1:2%:3
    for iter=pixels
        x_t=X(:,iter);
        %     clusts=randi(num_clusts);
        %     x_t=base(:,clusts);
        %     x_t=x_t+normrnd(0,.1,size(x_t));
        %     newX=[newX x_t];
        if mod(iter,100)==0
            disp(iter)
            %imagesc(W)
            %drawnow
        end
        y_old=Inf;
        count=0;
        while norm(y_old-y_t)>tolerance && count<500
            y_old=y_t;
            for i=1:numel(y_t)
                if ismember(i,active_set)
                    y_t(i)=max(W(i,:)*x_t - M(i,:)*y_t,0);
                end
            end
            count=count+1;
        end
        
        for i=1:numel(y_t)
            if ~ismember(i,active_set)
                delta_y=norm(x_t,2)^2-sum(y_t.^2)-y_t(i).^2;                
                if delta_y^2<=lambda || delta_y<0                   
                    y_t(i)=0;
                else
                    disp(delta_y^2)
                    y_t(i)=sqrt(delta_y);
                    active_set=[active_set i];
                    disp('add component')
                end
            end
        end
        
        for i=1:m
            y_t_hat(i)=y_t_hat(i)+y_t(i).^2;
            if ismember(i, active_set)
                for j=1:n
                    if w_nonneg
                        W(i,j) = max(0,W(i,j)+y_t(i)*(x_t(j)-W(i,j)*y_t(i))/y_t_hat(i));
                    else
                        W(i,j)= W(i,j)+y_t(i)*(x_t(j)-W(i,j)*y_t(i))/y_t_hat(i);
                    end
                end
                
                for j=1:m
                    if i~=j
                        M(i,j)= M(i,j)+y_t(i)*(y_t(i)-M(i,j)*y_t(i))/y_t_hat(i);
                    end
                end
            end
        end
        [~,cl_found]=max(y_t);
        detected_clusters=[detected_clusters cl_found];
        %disp([y_t; cl_found; clusts(iter)]')
        y_ts=[y_ts y_t];
        
    end
end
%
close all
for kk=1:size(W,1)
    subplot(2,1,1)
    mask=reshape(W(kk,:),[d1 d2]);    
    imagesc(mask)
    axis image    
    subplot(2,1,2)
    hold off
    plot(y_ts(kk,end-T:end))  
    hold all
    thr=mean(mask(:))+1*std(mask(:));
    plot(nansum(X.*repmat(reshape(mask',[d1*d2,1]),[1,T]),1))
    
    pause
end
%%
reord=[];
for p=pixels
    reord(pixels(p))=p;
end
%%
clusters=detected_clusters(end-d1*d2-1:end);
[~,bb,cc]=find(pixels);
classes=unique(clusters);
colormap(hot)

for c = classes
    mask=zeros(size(clusters));
    mask(clusters==c)=1;
    imagesc(reshape(mask(reord),[d1 d2]))
    pause
end
%%

%%
nclust=3;
a=kmeans(x',nclust,'distance','correlation')
imagesc(reshape(a,[d1 d2]))
%%
for kk=1:nclust
    imagesc(reshape(a.*(a==kk),[d1 d2]))
    pause
end