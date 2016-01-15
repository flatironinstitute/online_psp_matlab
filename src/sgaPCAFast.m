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
       U=U+gamma*(y'*x - tril(y'*y)*U);
%        for i=1:q
%              DU(i,:)=y(i)*x-y(i)*(y(1:i)*U(1:i,:));
%        end

    case 'SGA'
        U=U + gamma*(y'*x - 2*tril(y'*y)*U + diag(diag(y'*y,0))*U);
%         for i=1:q
%             U(i,:)=U(i,:)+gamma*(y(i)*x - 2*y(i)*(y(1:i)*U(1:i,:)) + y(i)*y(i)*U(i,:));
%         end
        
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
