function [eigval, eigvect] = CCovIncrPCA_fast(lambda, U, x, n, options)
% [eigval, eigvect] = incrPCA_fast(lambda, U, x, options)
% implementaiton of incermental PCA from onlinePCA R package
% parameters
% lambda: eigenvalues
% U: eigenvectors
% x: current sample
% options.n: sample size before observing x
% options.f: forgetting factor: a number in (0,1).
% options.q: number of components to find
% options.tol: tolerance
% output
% eigval:updated eigenvalues
% eigvect=iupdated eigenvectors

n=options.l;
q=options.q;

if isfield(options,'center')
    x = x - center;
end

if isfield(options,'tol')
    tol=options.tol;
else
    tol=1e-7;
end    



k =length(lambda);
q =int32(q);

enable_checks=1;

if enable_checks
    if (size(U,2) ~= k)
        error('length(lambda) ~= ncol(U)')
    end
    
    if (size(U,1) ~= length(x))
        error('length(x) ~= nrow(U)')
    end
    
   
end

% 
% v(:,i) = ((n-1-l)/n) * v(:,i) + ((1+l)/n) * u(:,i) * u(:,i)' * (v(:,i)/norm(v(:,i)));      % Euclidean length 
% u(:,i+1) = u(:,1,i) - u(:,i)' * (v(:,i)/norm(v(:,i))) * (v(:,i)/norm(v(:,i)));             

x1=x;

for i=1:min(q,n)
    U(:,i)=(n-1-l)/n*U(:,i)+(1+l)/n*x1*(x1'*(U(:,i)/norm(U(:,i))));
    x=x-x'*(U(:,i)/norm(U,:,i))*(U(:,i)/norm(U,:,i));    
end

% [eigval,ordr]=sort(eigval,'descend');   
% eigvect=V(:,ordr);
% end


% if (q<k)
%     eigval=eigval(1:q);
%     eigvect = eigvect(:,1:q);
% end

eigvect=U;

