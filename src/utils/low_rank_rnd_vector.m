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
% method: 'brownian_motion': covariance matrix in the form c_ij=max(i,j)
%         'spiked_covariance_normalized':  principal eigenvalues have the profile rho linspace(1,lambda_q,q)
%         'spiked_covariance': principal eigenvalues have the profile rho+ gap + slope*[q-1:-1:0]'
%         'ORL': use the orl database
%         'MNIST': use the MNIST database
% options: if 'spiked_covariance': 
%                     options.rho: size of non principal eigenvalues (i.e. noise)
%                     options.gap: difference between rho and first eigenvalue
%                     options.slope: slope with which the eignvalue grows departing from rho+gap
%          if 'spiked_covariance_normalized': 
%                     options.rho: size of non principal eigenvalues (i.e. noise) 
%                     options.lambda_q: size of smallest principal eigenvalues  
% 
% Returns
% --------
% X: vector of generated samples ( n x d )
% eig_vect: eigenvectors used for generating the datasets (d x q)
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
if isequal(method,'spiked_covariance') || isequal(method,'spiked_covariance_normalized')
    if isequal(method,'spiked_covariance')
        slope=options.slope;
        rho=options.rho;
        gap=options.gap;
        
        sigm = sqrt(gap + slope*[q-1:-1:0]');
    else
        rho=options.rho;
        lambda_q=options.lambda_q;        
        sigm = sqrt(linspace(1,lambda_q,q))';
    end
    
    eig_vect = orth(normrnd(0,1,d,q)); %compute orthonormal basis


    %sigm(1:q)=sqrt([1 .9 .8 .7 .6]);
    X=zeros(n,d);
    %disp('created X')
    x = normrnd(0,1,q,n);
    %disp('created x')
    eta = normrnd(0,1,d,n);
    %disp('creating data samples....')
    for kk=1:n
        if mod(kk,100)==0
            %disp(kk)
        end
%         X(kk,:) = sum(eig_vect*diag(sigm.*x(:,kk)),2)+sqrt(rho)*eta(:,kk);
        X(kk,:)=sum(bsxfun(@times,eig_vect,(sigm.*x(:,kk))'),2)+sqrt(rho)*eta(:,kk);
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
elseif isequal(method,'ORL') || isequal(method,'MNIST') || isequal(method,'YALE') || isequal(method,'ATT')
    
    if isequal(method,'ORL')        
        load('ORL_32x32');
        x1=fea';                 
        %load('pca_ORL') % load data eignvalues and eigenvectors
    elseif isequal(method,'MNIST')  
        x1 = loadMNISTImages('train-images-idx3-ubyte');
        %load('pca_MNIST') % load data eignvalues and eigenvectors
    elseif isequal(method,'YALE')
        load('YaleB_32x32.mat')
        x1=fea'; 
    elseif isequal(method,'ATT')
        load('ATT_faces_112_92.mat')
        x1=fea';     
    end
    mu=mean(x1,2);
    x=bsxfun(@minus,x1,mu); 
    
%     [eig_vect,~,eig_val]=pca(x','NumComponents',q);
%     eig_val=eig_val(1:q);
%     eig_vect=eig_vect(:,1:q);
    eig_vect=[];
    eig_val=[];  
    if n<=size(x,2)
        X=x(:,randperm(n))';
    else        
        X=x(:,randsample(1:size(x,2),n,1))';
    end
elseif isequal(method,'brownian_motion')
    
    x1 = normrnd(0,1/sqrt(d),n,d);
    
    X = cumsum(x1,2);    
    
     mu=mean(X,2);
     X=bsxfun(@minus,X,mu);  
%     
    if options.compute_eig
        
        if d^2<10e+8
            
            [I,J] = ind2sub([d,d],1:d^2);
            C=zeros(d);
            C(:)=min(I,J)/d;            
            [eig_vect,S]=eigs(C,q);
            eig_val=diag(S);
        
        else
            
            warning('eigenvalues not computed since covariace matrix contains more than 10e+7 elements')
            eig_vect=[];
            eig_val=[];
    
        end
        
    else
        eig_vect=[];
        eig_val=[];    
    end
    
else
    
    error('Undefined Method')

end
