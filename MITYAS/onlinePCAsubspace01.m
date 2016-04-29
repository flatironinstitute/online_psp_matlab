clear all
addpath('../MNIST')
% close all
q=64;
fea=0; % if 1 uses orl otherwise MNIST
if fea
    load('../MNIST/ORL_32x32');
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
    x1 = loadMNISTImages('../MNIST/train-images-idx3-ubyte');
    x1=x1(:,1:5000);
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

%%
X=x;
T=size(X,2);
in_dims=size(X,1);
out_dims=30;
w=zeros(out_dims,in_dims);
m=zeros(out_dims,out_dims);
sumy=zeros(out_dims,1);
Ysq=zeros(out_dims,1);
y=zeros(out_dims,T);
act=0;
eta=0.1;
cc=zeros(1,T);
tol=1e-5;
max_res=3e-1;
y_old=[];
lambda=0;
%%
tic
for i=1:T
    if mod(i,100)==0
        disp([i,act])
    end
    y_old(1:act,1)=-ones(act,1);

    tmp=sum((y(1:act,i)-y_old(1:act,1)).^2);
    
    y(1:act,i)=(eye(act)+m(1:act,1:act))\(w(1:act,:)*X(:,i));
    
    res=(sum(X(:,i).^2)-sum(y(1:act,i).^2))/sum(X(:,i).^2);
    if res>max_res
        if act<out_dims
            act=act+1;
            y(act,i)=sqrt(res*sum(X(:,i).^2));
            y_fast(1:act,i)=sqrt(res*sum(X(:,i).^2));
        else
            %disp('*')
        end
    end

    
    Ysq(1:act) = Ysq(1:act) + y(1:act,i).^2;

    Y_tmp=zeros(q,1);
    Y_tmp_sq=zeros(q,1);
    Y_tmp(1:act)=(y(1:act,i)./Ysq(1:act));
    Y_tmp_sq(1:act)=y(1:act,i).^2./Ysq(1:act);
    
    w(1:act,:) = w(1:act,:) +  bsxfun(@times,Y_tmp(1:act),X(:,i)') - bsxfun(@times,w(1:act,:),Y_tmp_sq(1:act));
    m(1:act,1:act) =m(1:act,1:act) + (1+lambda)*bsxfun(@times,Y_tmp(1:act), y(1:act,i)') - bsxfun(@times,m(1:act,1:act),Y_tmp_sq(1:act));
    m=m-diag(diag(m));
end
toc
% %%
% tic
% for i=1:T
%     if mod(i,100)==0
%         disp([i,act])
%     end
%     y_old(1:act,1)=-ones(act,1);
%     cc(i)=0; %?
%     tmp=sum((y(1:act,i)-y_old(1:act,1)).^2);
%     
%     while cc(i)<100 && tmp>tol
%         y_old=y(1:act,i);
%         y(1:act,i)=y_old(1:act,1)*(1-eta)+eta*(w(1:act,:)*X(:,i)-m(1:act,1:act)*y_old(1:act,1));
%         cc(i)=cc(i)+1;
%     end
%     res=(sum(X(:,i).^2)-sum(y(1:act,i).^2))/sum(X(:,i).^2);
%     if res>max_res
%         if act<out_dims
%             act=act+1;
%             y(act,i)=sqrt(res*sum(X(:,i).^2));
%         else
%             %disp('*')
%         end
%     end
%     sumy(1:act,1)=sumy(1:act,1)+y(1:act,i).^2;
%     w(1:act,1:in_dims)=w(1:act,1:in_dims)+(y(1:act,i)*X(:,i)'-diag(y(1:act,i).^2)*w(1:act,1:in_dims))./(sumy(1:act,1)*ones(1,in_dims));
%     m(1:act,1:act)=m(1:act,1:act)+(y(1:act,i)*y(1:act,i)'-diag(y(1:act,i).^2)*m(1:act,1:act))./(sumy(1:act,1)*ones(1,act));
%     m=m-diag(diag(m));
% end
% toc
%%
%figure,imagesc(w)
figure
f=(m+eye(out_dims))\w;
figure
for i=1:act
    subplot(10,20,i),imagesc(reshape(f(i,:),28,28))
end