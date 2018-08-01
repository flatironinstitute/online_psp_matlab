function [X, U, lam ] = standardize_data(X, U, lam)
% demean and 
X = bsxfun(@minus, X , mean(X,2));
scale_factor  = 1/mean(sqrt(sum(X.^2, 1)));
X = X * scale_factor;
lam = lam * (scale_factor^2);