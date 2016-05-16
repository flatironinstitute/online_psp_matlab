clear all
addpath('../MNIST')
% close all

x1 = loadMNISTImages('../MNIST/train-images-idx3-ubyte');
x1=x1(:,1:end);
d1=28;
d2=28;
x=x1;
%     mu=mean(x1,2);
%     x=bsxfun(@minus,x1,mu);

total_img=size(x,2);

%%
X=x;
T=size(X,2);
in_dims=size(X,1);
out_dims=30;
w=zeros(out_dims,in_dims);
m=zeros(out_dims,out_dims);
Ysq=zeros(out_dims,1);
y=zeros(out_dims,T);
act=0;
eta=0.1;
max_res=1e-1;
y_old=[];
lambda=0;
%
tic
for i=1:T
    disp(act)
    disp(['*' num2str(sum(Ysq>0))])
    if mod(i,100)==0
        disp([i,act])
    end
    
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

%%
%figure,imagesc(w)
figure
f=(m+eye(out_dims))\w;

    
for i=1:act
    subplot(5,6,i),imagesc(reshape(f(i,:),d1,d2))
end