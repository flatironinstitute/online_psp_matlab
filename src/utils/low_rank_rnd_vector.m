function [X,eig_vect,eig_val] = low_rank_rnd_vector(d,q,n,method,options)

% generate random samples corresponding to a particular eigenvalue
% configuration. 
% 
% 
% Parameter
% ---------
% d: input dimensionality 
% q: number of principal components
% n: number of samples to generate 
% 
% method: 'brownian_motion' or 'spiked_covariance'
% 
% options: if 'spiked_covariance': 
%                     options.rho: size of non principal eigenvalues (i.e. noise)
%                     options.gap: difference between rho and first eigenvalue
%                     options.slope: slope with which the eignvalue grows departing from rho+gap
% 
% 
% Returns
% --------
% X: vector of generated samples ( n x d )
% eig_vect: eigenvectors used for generating the datasets
% eig_val: corresponding original eigenvalues such that the stimated
% eigenvalues should be [eig_val+rho; rho*ones(d-q,1)] for 'spiked_covariance'
% 
% Example
% --------
%
% d=50;
% q=5;
% n=2000;
% 
% method='spiked_covariance';
% 
% options.rho=0.1;
% options.gap=0.4;
% options.slope=.1;
% 
% [X,eig_vect,eig_val] = low_rank_rnd_vector(d,q,n,method,options);
%
% [COEFF, SCORE, LATENT] = pca(X,'NumComponents',q);
% disp('PCA')
% plot(LATENT)
% compute_reconstruction_error(eig_vect,COEFF)
% 
% close all
% subplot(2,2,1)
% hold all
% plot([eig_val+options.rho; options.rho*ones(d-q,1)])
% plot(LATENT)
% subplot(2,2,2)
% imagesc(eig_vect*eig_vect')
% subplot(2,2,3)
% imagesc(COEFF*COEFF')
% subplot(2,2,4)
% imagesc(abs(COEFF*COEFF'-eig_vect*eig_vect'),[0 .01])

%%
if isequal(method,'spiked_covariance')
    slope=options.slope;
    rho=options.rho;
    gap=options.gap;

    eig_vect = orth(normrnd(0,1,d,q)); %compute orthonormal basis
    sigm = sqrt(gap + slope*[q-1:-1:0]');

    %sigm(1:q)=sqrt([1 .9 .8 .7 .6]);
    X=zeros(n,d);
    %disp('created X')
    x = normrnd(0,1,q,n);
    %disp('created x')
    eta = normrnd(0,1,d,n);
    disp('creating data samples....')
    for kk=1:n
        if mod(kk,100)==0
            disp(kk)
        end
        X(kk,:)=sum(eig_vect*diag(sigm.*x(:,kk)),2)+sqrt(rho)*eta(:,kk);
    end
    
    eig_val=sigm.^2;
    
    % C=cov(X');
    % disp('created C')
    % 
    % [vectors,lam] = eig(C,'vector');
    % disp('EIG')
    
    % %vectors=vectors(end:-1:end-4,:)';
    % vectors=vectors(:,end:-1:end-4);
    % [U,S,V] = svd(C,'econ');
    % disp('SVD')
    
else
    error('Undefined Method')
end
