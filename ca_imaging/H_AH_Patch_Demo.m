clear all
nam = 'M_FLUO.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=1452;					% user input: how many frames to read   (optional, default until the end)

mov = bigread2(nam,sframe,num2read);
if ~isa(mov,'double');    mov = double(mov);  end         % convert to double
mov=mov(10:end-10,10:end-10,:);
mov=mov(:,30:60,:);
[d1,d2,T]=size(mov);
M=reshape(mov,[d1*d2 T])';
% imagesc(Yr./max(Yr(:)))
M(isnan(M))=0;
%[wcoeff,score,latent,tsquared,explained]  = princomp(M);
M=bsxfun(@minus,M,quantile(M,.5,1));
x=M;
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
sigma = cov(x');% x * x' / size(x, 2);
k=40;
[U,S,V] = svds(sigma,k);
epsilon=0;
zca_whitesig =x'*U*diag(1./sqrt(diag(S) + epsilon))*U';
zca_whitesig =zca_whitesig';
%%
tic
sigma = cov(x');% x * x' / size(x, 2);
k=100;
[U,S,V] = svds(sigma,k);

xRot = U' * x;          % rotated version of the data.

xTilde = U(:,1:k)' * x; % reduced dimension representation of the data,
% where k is the number of eigenvectors to keep
epsilon=0;
zca_whitesig = U * diag(1./sqrt(diag(S) + epsilon)) * U' * x;
%zca_whitesig = diag(1./sqrt(diag(S) + epsilon)) * U' * x;


zca_whitesig = U * U' * x;

toc
rec_error=norm(x-U*U'*x,2)/norm(x,2)
%%
close all
X=x';
X=bsxfun(@times,x',1./(std(x',[],2)));
%lambda=500000; % for patch
% X=(x')./(max(x(:))-min(x(:)));
%lambda=1.5;
%X=X/norm(X,'fro');
%X=zca_whitesig';
lambda=70e+5;
% clc
tic
T=size(X,2);
m=16; %max_components
n=size(X,1);

w_nonneg=1; % impose nonnegativity on W
tolerance=1e-9;

active_set=[];
W=zeros(m,n);
M=zeros(m,m);
y_t=zeros(m,1);
y_t_hat=zeros(m,1);
detected_clusters=[];
newX=[];
y_ts=[];
pixels=randperm(T); % if you want to present frames randomly
% pixels=1:T;
%
for nnn=1:5
    for iter=pixels
        x_t=X(:,iter);        
        if mod(iter,100)==0
            disp(iter)
%             imagesc(W)
%             drawnow
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
        y_t_hat=y_t_hat+y_t.^2;
        
%         if ~isempty(active_set)
%             W(active_set,:) = max(0,W(active_set,:)+y_t(active_set)*(x_t'-y_t(active_set)'*W(active_set,:))./repmat(y_t_hat(active_set),[1,n]));
%         end
        for i=1:m
            %y_t_hat(i)=y_t_hat(i)+y_t(i).^2;
             if ismember(i, active_set)
                if w_nonneg
                    W(i,:) = max(0,W(i,:)+y_t(i)*(x_t'-W(i,:)*y_t(i))/y_t_hat(i));
                else
                    W(i,:) = W(i,:)+y_t(i)*(x_t'-W(i,:)*y_t(i))/y_t_hat(i);
                end
                
                idx=setdiff(1:m,i);
                M(i,idx)= M(i,idx) +y_t(i)*(y_t(i)-M(i,idx)*y_t(i))/y_t_hat(i);
%                 for j=1:m
%                     if i~=j
%                         M(i,j)= M(i,j)+y_t(i)*(y_t(i)-M(i,j)*y_t(i))/y_t_hat(i);
%                     end
%                 end

              
            end
        end
%         
        if(mod(iter,50)==0)
            W1=reshape(W,[],d1,d2);        
            for lll=1:size(W1,1)
                if sum(sum(W1(lll,:,:)))>0
                    subplot(4,4,lll)
                    imagesc(squeeze(W1(lll,:,:)))
                end
            end
        drawnow
        end
        
        [~,cl_found]=max(y_t);
        detected_clusters=[detected_clusters cl_found];
        %disp([y_t; cl_found; clusts(iter)]')
        y_ts=[y_ts y_t];
        
    end
end
%%
figure
n_comp=12;
[A,B]=nnmf(x',n_comp);
for tt=1:n_comp
    subplot(4,3,tt)
    imagesc(reshape(A(:,tt),d1,d2))
end
%%
close all
for kk=1:size(W,1)
    subplot(2,1,1)
    mask=reshape(W(kk,:),[d1 d2]);    
    imagesc(mask)
    axis image    
    subplot(2,1 ,2)
    hold off
    plot(y_ts(kk,end-T:end))  
%     plot(y_ts(kk,:))  
 

    hold all
    thr=mean(mask(:))+1*std(mask(:));
    plot(nansum(X.*repmat(reshape(mask',[d1*d2,1]),[1,T]),1))
    
    pause
end