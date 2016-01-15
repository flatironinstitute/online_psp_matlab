function [eigval, eigvect] = incrPCA_fast(lambda, U, x, options)
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

n=options.n;
q=options.q;

if isfield(options,'center')
    x = x - center;
end

if isfield(options,'tol')
    tol=options.tol;
else
    tol=1e-7;
end    

if isfield(options,'f')
    f=options.f;
else
    f=1/n;
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

x=x';
lambda = (1-f) * lambda;
x = sqrt(f) * x;
xhat =U' * x; 
x = x - U * xhat; % project to the orthogonal space
normx = sqrt(sum(x.^2)); % how largeis orthogonal space?
if (normx >= tol)
    k = k+int32(1); % increase vector by 1 ans set to 0
    lambda(k) = 0; 
    xhat(k) = normx;
    U = [U,x/normx];
    %dim(U) <- c(length(x),k)
end

% if 0
%     [V,D,~] = svd(diag(lambda)+xhat*xhat');
%     eigval=diag(D); 
%     eigvect=V;
% else 
[V,D,~] = eig(diag(lambda)+xhat*xhat');
eigval=diag(D);
[eigval,ordr]=sort(eigval,'descend');   
eigvect=V(:,ordr);
% end
    
if (q<k)
    eigval=eigval(1:q);
    eigvect = eigvect(:,1:q);
end

eigvect=U*eigvect;
