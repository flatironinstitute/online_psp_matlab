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
%%
for i=1:T
    [i,act]
        y_old(1:act,1)=-ones(act,1);
        cc(i)=0;
        %sum(y_old(1:act,1).^2)
    
        if act>0
            while cc(i)<1000 && sum((y(1:act,i)-y_old(1:act,1)).^2)>(1e-8)*sum(y_old(1:act,1).^2)
          %  while  (cc(i)<2000 && max(abs(y(1:act,i)-y_old(1:act,1)))>1e-6)
        
       y_old(1:act,1)=y(1:act,i);
       y(1:act,i)=max(y_old(1:act,1)*(1-eta)+eta*(w(1:act,:)*X(:,i)-m(1:act,1:act)*y_old(1:act,1))-th(1:act,i),0);
       cc(i)=cc(i)+1;
            end
        end
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
%%
figure,imagesc(w)
f=(m+eye(out_dims))\w;
figure
for i=1:act
    subplot(10,20,i),imagesc(reshape(f0(i,:)*U(:,:)',28,28))
%subplot(7,13,i),imagesc(reshape(f(i,:)*U(:,1:100)',12,12))
end