% number of examples
n=size(X,2);
% input dimensions
in_dims=size(X,1);
% output dimensions
out_dims=1000;
% offdiagonal cost term promoting orthogonalization
gamma=3;
% cost for adding a new output channel
delta=0.9;
%delta=0.98;
% neuronal threshold defining locality in data space
theta=0.;
% initialization of synaptic weights
w=zeros(out_dims,in_dims);
m=zeros(out_dims,out_dims);
% initialization of cumulative activity, nonzero value seems to reduce the
% impact of the first activating data point
sumy=1*ones(out_dims,1);
y=zeros(out_dims,n);
% number of active output channels
act=0;
% single pass over data
for i=1:n
    [ i act ]
    % visualization of the feature vectors
    if rem(i,100)==0
        f=(m+eye(out_dims))\w;
        close,figure
        for ii=1:act
            hold on
            %plot(X(1,1:i),X(2,1:i),'r.')
            subplot(30,40,ii),imagesc(reshape(f(ii,:),28,28))
            %plot(f(ii,1),f(ii,2),'o')
            drawnow
        end
    end
    % initialize output
    y_old(1:act,1)=-ones(act,1);
    % coordinate descent iterations
    for k=1:1000
        if sum((y(1:act,i)-y_old(1:act,1)).^2)/sum((y_old(1:act,1)).^2)>1e-6
            y_old=y(:,i);
            %y(1:act,i)=w(1:act,:)*X(:,i)-m(1:act,1:act)*y_old(1:act,1);
            % coordinate descent loop
            for t=1:act
                % defining a cubic equation
                c=[4 0 2*(delta+2*(y_old'*y_old-y_old(t,1)^2)-2*X(:,i)'*X(:,i))+4*sumy(t) 4*sumy(t)*(m(t,1:act)*y_old(1:act,1)-w(t,:)*X(:,i))];
                r=roots(c);
                % finding a real root
                nn=find(abs(imag(r))<1e-3);
                %y(t,i)=max(max(r(nn)),0);
                % finding greatest real root
                y(t,i)=max(r(nn));
                % hard-thresholding output
                y(t,i)=y(t,i)*max(y(t,i)-theta,0)/(y(t,i)-theta);
            end
            %y(1:act,i)=(y(1:act,i)+abs(y(1:act,i)))/2;
        end
    end
    %activating another output channel
    if (sum(y(1:act,i).^2)-sum(X(:,i).^2))/sum(X(:,i).^2)<-delta
        act=act+1;
        %y(act,i)=norm(X(:,i)-w(1:act,:)'*y(1:act,i));
        y(:,i)=zeros(out_dims,1);y(act,i)=norm(X(:,i));
        %y(act,i)=sqrt(sum(X(:,i).^2)-delta/2-sum(y(1:act,i).^2));
    end
    % updating cumulative output activity
    sumy(1:act,1)=sumy(1:act,1)+y(1:act,i).^2;
    % updating synaptic weights
    w(1:act,1:in_dims)=w(1:act,1:in_dims)+(y(1:act,i)*X(:,i)'-diag(y(1:act,i).^2)*w(1:act,1:in_dims))./(sumy(1:act,1)*ones(1,in_dims));
    m(1:act,1:act)=m(1:act,1:act)+((1+gamma)*y(1:act,i)*y(1:act,i)'-diag(y(1:act,i).^2)*m(1:act,1:act))./(sumy(1:act,1)*ones(1,act));
    m=m-diag(diag(m));
end
%figure,imagesc(w)
% f=(m+eye(10))\w;
% figure
% for i=1:10
%     subplot(2,5,i),imagesc(reshape(f(i,:),28,28))
% end