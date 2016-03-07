function [eigval, eigvect] = sgaPCAFast(lambda, U, x,options)
% [eigval, eigvect,y] = sgaPCAFast(lambda, U, x,options)
% implementaiton of incermental PCA from onlinePCA R package
% inputs:
% lambda: eigenvalues
% U: eigenvectors ()
% x: current sample
% options.gamma: learning rate 
% options.q: number of components to find
% options.method:'SGA' (stochastic gradient ascend) or GHA (General Hebbian Algorithm)'
% options.do_sort: boolean, whether to sort the components or not
% outputs
% eigval: updated eingenvalues
% eigvect: updated eigenvectors

gamma=options.gamma;
q=options.q;
do_sort=options.do_sort;
method=options.method;
d=size(U,1);
k = size(U,2);
enable_checks=0;

if enable_checks
    
    if length(x) ~= d
        error('Dimensions of X not matching d')
    end
    
    if exist('lambda')
        if length(lambda) ~= k
            error('dimensions of lambda not matching k')
        end
    end
   
end

y=(x*U);
U=U';
%DU=zeros(size(U));
% apply update rule on eigenvectors. 
switch method
    case 'GHA'
       % DU=y'*x + U*tril(y*y');
%        U=U+gamma*(y'*x - tril(y'*y)*U);
         error('work in progress')

    case 'SGA'
       phi=U*x';
       res1=bsxfun(@times,phi,U);
       cum_res=cumsum(res1(1:end-1,:));
       res=[zeros(1,d); cum_res];
       U = U + gamma*bsxfun(@times,bsxfun(@plus,-res1-2*res,x),phi);
%         U2=U + gamma*(y'*x - 2*tril(y'*y)*U + sparse(diag(diag(y'*y,0)))*U);       

        
    otherwise
        error 'unknown method'
end

U=U';
y=y';
y=y(1:q);

lambda=(1-gamma)*lambda+gamma*y.*y;    
if do_sort
    [~,ix] = sort(lambda,1, 'descend');
    if ~isequal(ix,1:k)
        lambda = lambda(ix);
        U = U(:,ix);
    end       
end
if q < k 
    lambda=lambda(1:q);
end

eigval=lambda;
eigvect=U;
