function [values,vectors,iter_final]=run_incrPCA(x,q,n_init_PCA,values,vectors,iter_so_far,options_algorithm)

if ~exist('options_algorithm')
    options_algorithm=struct();
    options_algorithm.pca_algorithm='IPCA';    
    options_algorithm.tol=1e-7;
end

options_algorithm.q=q;
[d,T]=size(x);
scramble=randperm(T);
x=x(:,scramble);

    
if n_init_PCA>0
    if n_init_PCA<q
        error('n_init_PCA cannot be smaller than q')
    end
    if ~isempty(vectors) ||  ~isempty(values) ||  ~isempty(iter_so_far)
        error('you cannot have n_init_PCA >0 and pass inputs vectors, values , iter_so_far')
    end

    n0=n_init_PCA;
    [eigvect_init,~,eigval_init]=pca(x(:,1:(n0+1))','NumComponents',q);     
    vectors=eigvect_init; 
    values=eigval_init(1:q); 
else
    if isempty(values) || isempty(vectors) 
        n0=1;
        values=randn(q,1)/sqrt(d);
        vectors=randn(d,q)/sqrt(d);
        vectors(:,1)=x(:,1)/norm(x(:,1));
        values(1)=max(x(:,1));
        iter_so_far=0;
        disp('Initializing eigenvalues and eigenvectors randomly')    
    end
end


if isempty(iter_so_far)
   iter_so_far=0; 
end


for iter=(n_init_PCA+1):(T)    
    if mod(iter,1000)==0
        disp(iter_so_far+iter-1)
    end
    options_algorithm.n=max(1,iter_so_far+iter-1);
    
    [values, vectors] = incrPCA_fast(values, vectors, x(:,iter)', options_algorithm);        
end

iter_final=iter;


