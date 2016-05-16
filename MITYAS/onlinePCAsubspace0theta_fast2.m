clear all
format compact
addpath('../MNIST')
% close all

x1 = loadMNISTImages('../MNIST/train-images-idx3-ubyte');
x1=x1(:,1:end);
d1=28;
d2=28;
X=x1;
mu=mean(x1,2);
x=bsxfun(@minus,x1,mu);
mu=mean(x,2);
x=bsxfun(@minus,x,mu);
total_img=size(X,2);
%%
T=size(X,2);
in_dims=size(X,1);
out_dims=1000;
w=zeros(out_dims,in_dims);
m=zeros(out_dims,out_dims);
sumy=zeros(out_dims,1);
y=zeros(out_dims,T);
act=0;
eta=0.1;
cc=zeros(1,T);
theta=0;
th=zeros(out_dims,T);
%
tic
for i=1:2000%T
    if mod(i,100)==0
        disp([i,act])
    end
%     
    y_old(1:act,1)=-ones(act,1);
    cc(i)=0;
    %sum(y_old(1:act,1).^2)
    st=toc;
%     disp('*')
    if act>0
% %         while cc(i)<1000 && sum((y(1:act,i)-y_old(1:act,1)).^2)>(1e-8)*sum(y_old(1:act,1).^2)
%         while  (cc(i)<10000 && max(abs(y(1:act,i)-y_old(1:act,1)))>1e-5)            
%             y_old(1:act,1)=y(1:act,i);
%             y(1:act,i)=max(y_old(1:act,1)*(1-eta)+eta*(w(1:act,:)*X(:,i)-m(1:act,1:act)*y_old(1:act,1))-th(1:act,i),0);
%             cc(i)=cc(i)+1;
%         end
        y(1:act,i) = NSM_neuraldynMat(w(1:act,:)*X(:,i),m(1:act,1:act),1e-5,rand(size(y_old(1:act,1))));
        
    end
%     aa=(norm(y(1:act,i)-w(1:act,:)*X(:,i)+m(1:act,1:act)*y(1:act,i)));
    y_scd=sum(y(1:act,i).^2);
%     a1=toc-st;
%     if mod(i,10)==0 && i>100
        
%         n1=norm(X(:,i-100:i)'*X(:,i-100:i)-y(1:act,i-100:i)'*y(1:act,i-100:i),'fro');
%         disp(n1)
%     end      
%         %     disp(n1)
%     st=toc;
%        % disp(norm(w(1:act,:)*X(:,i)-(eye(act)+m(1:act,1:act))*y(1:act,i),'fro')/norm(y(1:act,i),'fro'))
      [y(1:act,i),resnorm] = lsqnonneg((eye(act)+m(1:act,1:act)),w(1:act,:)*X(:,i));
      disp(max(y(1:act,i)))
%       disp(norm(y(1:act,i)-w(1:act,:)*X(:,i)+m(1:act,1:act)*y(1:act,i))/aa)
%     a2=toc-st;
% %     
%     if mod(i,10)==0 && i>100
% %         disp(a2/a1)
%         n2=norm(X(:,i-100:i)'*X(:,i-100:i)-y(1:act,i-100:i)'*y(1:act,i-100:i),'fro');
% 
% %         n2=norm(X(:,1:i)'*X(:,1:i)-y(1:act,1:i)'*y(1:act,1:i),'fro')/norm(y(1:act,1:i)'*y(1:act,1:i),'fro');
%         %     disp(norm(w(1:act,:)*X(:,i)-(eye(act)+m(1:act,1:act))*y(1:act,i),'fro')/norm(y(1:act,i),'fro'))
%         disp(n2)
%     end
%      disp(y_scd/sum(y(1:act,i).^2))
    
    res=(sum(X(:,i).^2)-sum(y(1:act,i).^2))/sum(X(:,i).^2);
    if res>1e-2
        if act<out_dims
            act=act+1;
            y(act,i)=sqrt(res*sum(X(:,i).^2));
        else
            'increase out_dims'
        end
    end
    sumy(1:act,1)=sumy(1:act,1)+y(1:act,i).^2;
    w(1:act,1:in_dims)=w(1:act,1:in_dims)+(y(1:act,i)*X(:,i)'-diag(y(1:act,i).^2)*w(1:act,1:in_dims))./(sumy(1:act,1)*ones(1,in_dims));
    m(1:act,1:act)=m(1:act,1:act)+(y(1:act,i)*y(1:act,i)'-diag(y(1:act,i).^2)*m(1:act,1:act))./(sumy(1:act,1)*ones(1,act));
    m=m-diag(diag(m));
    th(1:act,i+1)=th(1:act,i)+(theta*y(1:act,i)-diag(y(1:act,i).^2)*th(1:act,i))./(sumy(1:act,1));
end
toc
%%
figure,imagesc(w)
f=(m+eye(out_dims))\w;
figure
for i=1:act
    subplot(10,20,i),imagesc(reshape(f0(i,:)*U(:,:)',28,28))
    %subplot(7,13,i),imagesc(reshape(f(i,:)*U(:,1:100)',12,12))
end