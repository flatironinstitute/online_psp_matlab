function [M,W,Ysq,Ys]=run_H_AH_PCA(x,q,n_init_PCA,W,M,Ysq,options_algorithm)

if ~exist('options_algorithm')
    options_algorithm=struct();
    options_algorithm.pca_algorithm='H_AH_NN_PCA';
    options_algorithm.update_method='ls';
    options_algorithm.tol=1e-5;
end



options_algorithm.q=q;
[d,T]=size(x);
% scramble=randperm(T);
% x=x(:,scramble);

if n_init_PCA>0
    if n_init_PCA<q
        error('n_init_PCA cannot be smaller than q')
    end
    if ~isempty(W) || ~isempty(M) || ~isempty(Ysq)
        error('you cannot have n_init_PCA >0 and pass inputs W,M,Ysq')
    end

    n0=n_init_PCA;
    [eigvect_init,~,eigval_init]=pca(x(:,1:(n0+1))','NumComponents',q);
    vectors=eigvect_init;
    M = zeros(q,q);
    W = vectors';
    Ysq=sum((W*x(:,1:n0)).^2,2);
else
    if isempty(W) && isempty(M) && isempty(Ysq)

        W=randn(q,d)/sqrt(d); 
        W(1,:)=x(:,1)/norm(x(:,1));
        M=zeros(q,q);
        M=M-diag(diag(M));
        %Ysq=10*ones(size(W,1),1);
        Ysq=2*pi*norm(x(:,1),'fro')^2/d*ones(size(W,1),1);
%         Ysq=2*pi*norm(x(:,1),'fro')^2/sqrt(d)*ones(size(W,1),1);

    end
end

Y=rand(q,1)*0;

% for iter=1:n_iter
%     disp(['** Iteration:' num2str(iter)])
%     scramble=randperm(T);
%     for kk=1:T
%         if mod(kk,1000)==0
%             disp(kk)
%         end
%         options_algorithm.gamma=1./Ysq;
%         [M,W,Y]=H_AH_NN_PCAFast(M,W,Y,x(:,scramble(kk))',options_algorithm);
%         Ysq = Ysq + Y.^2;
%         Y_out(:,kk)=Y;
%     end
% end
Ys = [];
for iter=(n_init_PCA+1):(T)    
    if mod(iter,100)==0
        disp(iter)
    end
%     options_algorithm.gamma=1./Ysq;
    [M,W,Y,Ysq]=H_AH_NN_PCAFast(M,W,Y,Ysq,x(:,iter)',options_algorithm);
     Ys = [Ys Y];     
     
end
disp(iter)


