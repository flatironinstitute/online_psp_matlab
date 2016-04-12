function [errors_real,errors_batch_pca,errors_online,times_,ff_name]=Online_PCA_simulations(folder_exp,options_simulations,options_generator,options_algorithm)
%%
method_random=options_generator.method;

q=options_generator.q;
d=options_generator.d;
n=options_generator.n;

outer_iter=options_simulations.outer_iter;
n0=options_simulations.n0;
niter=options_simulations.niter;
nstep_skip_EIGV_errors=options_simulations.nstep_skip_EIGV_errors;
compute_error_batch=options_simulations.compute_error_batch;
compute_error_real=options_simulations.compute_error_real;
compute_error_online=options_simulations.compute_error_online;
normalize_input=options_simulations.normalize_input;
normalization_scale=options_simulations.normalization_scale;
initialize_PCA=options_simulations.initialize_PCA;
orthonormalize_vectors=options_simulations.orthonormalize_vectors;
init_weight_mult=options_simulations.init_weight_mult;

pca_algorithm=options_algorithm.pca_algorithm;
set(0, 'DefaulttextInterpreter', 'none')
if d>q
    hold all;
    errors_real=nan*zeros(n*outer_iter,niter);
    errors_reconstr=nan*zeros(1,niter);
    errors_batch_pca=nan*zeros(n*outer_iter,niter);
    errors_online=nan*zeros(n*outer_iter,niter);
    errors_ortho=nan*zeros(n*outer_iter,niter);
    errors_decorr=nan*zeros(n*outer_iter,niter);
    times_=nan*zeros(n*outer_iter,niter);
    for ll=1:niter
        disp(ll)
        Cy=[];
        %generate random samples
        if isequal(method_random,'brownian_motion') && ll>1 % eigenvalues are only computed once since the covariance matrix is always the same
            options_generator.compute_eig=0;
            [x,~,~] = low_rank_rnd_vector(d,q,n,method_random,options_generator);
        else
            if compute_error_online || compute_error_batch || compute_error_real
                options_generator.compute_eig=1;
            else
                options_generator.compute_eig=0;
            end
            [x,eig_vect_real,eig_val_real] = low_rank_rnd_vector(d,q,n,method_random,options_generator);
            mu=mean(x,2);
            x=bsxfun(@minus,x,mu);
            if normalize_input                
               x=normalization_scale*x; 
%                x=x/range(x(:));
               
%                  x=2*(x-min(x(:)))/(max(x(:))-min(x(:)))-1;   
%                  mu=mean(x,2);
%                  x=bsxfun(@minus,x,mu);
%                  x=x*sqrt(q);
               % x=x/quantile(abs(x(:)),.95);
%                 x=x/3/std(x(:));                
%                 mu=mean(x,2);
%                 x=bsxfun(@minus,x,mu);
%                 x=x/max(abs(x(:)));
            end
            if compute_error_batch
                [eig_vect_batch_pca,~,eig_val_batch_pca]=pca(x,'NumComponents',q);
                eig_val_batch_pca=eig_val_batch_pca(1:q);
            end
        end
        
        if initialize_PCA
            n0=max(n0,q);
            [eigvect_init,~,eigval_init]=pca(x(1:(n0+1),:),'NumComponents',q);            
        else
            if isequal('H_AH_NN_PCA',pca_algorithm) || isequal('GHA',pca_algorithm) || isequal('SGA',pca_algorithm)
                n0=1;
                eigval_init=randn(q,1)/sqrt(q);
                eigvect_init=randn(d,q)/sqrt(d);
                eigvect_init(:,1)=x(1,:)/norm(x(1,:));
%                 eigval_init(1)=max(x(1,:));
            else
                n0=1;
                eigval_init=randn(q,1)/sqrt(q);
                eigvect_init=randn(d,q)/sqrt(d);
                eigvect_init(:,1)=x(1,:)/norm(x(1,:));
                %eigval_init(1)=max(x(1,:));
            end
            
        end
        
        values=eigval_init(1:q);
        vectors=eigvect_init;
        
        if isequal('H_AH_NN_PCA',pca_algorithm)
            M = zeros(q,q);
            W = vectors';
            Y = zeros(q,1);
            %Ysq=max(10*ones(size(W,1),1),sum((W*x(1:n0,:)').^2,2));
            Y_=(W*x(1:n0,:)');
            Cy=Y_*Y_';
            if initialize_PCA
                Ysq=sum((W*x(1:n0,:)').^2,2);
            else
                Ysq=init_weight_mult*norm(x(1,:),'fro')^2/d*ones(size(W,1),1);
            end
        elseif isequal('GHA',pca_algorithm) || isequal('SGA',pca_algorithm)
            learning_rate=init_weight_mult;
        end
        
        
        if compute_error_batch || compute_error_real || compute_error_online
            if orthonormalize_vectors
                vectors_err = orth(vectors);
            end
        end
        
        if compute_error_online
            errors_online(n0,ll)=compute_projection_error(eig_vect_batch_pca,vectors_err);
        end
        
        if compute_error_batch
            errors_batch_pca(n0,ll)=compute_projection_error(eig_vect_batch_pca,vectors_err);
        end

        if compute_error_real
            errors_real(n0,ll)=compute_projection_error(eig_vect_real,vectors_err);
        end
            
        
        
        
        
        tic
        for outit=1:outer_iter
            for i = (n0+1):n
                idx=(outit-1)*n + i;
                
                switch pca_algorithm
                    case 'SGA'
                        options_algorithm.gamma=learning_rate/idx;
                        [values, vectors] = sgaPCAFast(values, vectors, x(i,:),options_algorithm);
                    case 'GHA'
                        options_algorithm.gamma=learning_rate/idx;
                        [values, vectors] = sgaPCAFast(values, vectors, x(i,:),options_algorithm);
                    case 'IPCA'
                        options_algorithm.n=idx-1;
                        [values, vectors] = incrPCA_fast(values, vectors,  x(i,:), options_algorithm);
                    case 'H_AH_NN_PCA'                        
                        [M,W,Y,Ysq] = H_AH_NN_PCAFast(M,W,Y,Ysq,x(i,:),options_algorithm);
                end
                
                if (numel(find(isnan(vectors)))==0) && (compute_error_online + compute_error_batch + compute_error_real) && (mod(idx,nstep_skip_EIGV_errors) ==  0 || idx==n*outer_iter || i == round(n*outer_iter/2))
                    if compute_error_online                    
                        [eig_vect_online,~,eigval_online]=pca(x(1:idx,:),'NumComponents',q);
                        eigval_online=eigval_online(1:q);                        
%                     else
%                         eig_vect_online=NaN;
                    end
                    

                    
                    if isequal('H_AH_NN_PCA',pca_algorithm) 
                        F=(pinv(diag(ones(q,1))+M(1:q,1:q))*W(1:q,:))';
%                         errors_ortho(idx,ll) = (norm(F'*F-eye(q),'fro')/norm(F*F','fro'));                        
                        Cy=Cy+Y*Y'; 
                        errors_decorr(idx,ll)=10*log10(norm((Cy-eye(q))/i,'fro')^2);
                        if orthonormalize_vectors
                            vectors_err = orth(F);
                        else
                            vectors_err = F;
                        end                        
                    elseif isequal('SGA',pca_algorithm) || isequal('GHA',pca_algorithm)  
                        
                        if orthonormalize_vectors
                            vectors_err = orth(vectors);
                        else
                            vectors_err=vectors;
                        end
                                                
                    else
                        vectors_err=vectors;
                    end
                    
                    if compute_error_real
                        errors_real(idx,ll)=compute_projection_error(eig_vect_real,vectors_err);
                    end
                    if compute_error_batch
                        errors_batch_pca(idx,ll)=compute_projection_error(eig_vect_batch_pca,vectors_err);
                    end
                    if compute_error_online
                        errors_online(idx,ll)=compute_projection_error(eig_vect_online,vectors_err);
                    end

                else
                    times_(idx,ll)=toc;
                end

            end
        end
        errors_reconstr(ll)=norm(x'-vectors_err*vectors_err'*x','fro')^2/norm(x,'fro')^2;

    end
    hold on    
    colrand=rand(1,3);
    if sum(colrand)>2.5
        colrand=[.5 .5 .5];
    end
    
    
    drawnow
    ff_name=['n' num2str(n) '_d' num2str(d) '_q' num2str(q) '_iw' num2str(init_weight_mult) '_ns' num2str(normalization_scale)]; 
    if isfield(options_algorithm,'lambda')
           ff_name=[ff_name '_lambda' num2str(options_algorithm.lambda)];  
    end
    if isfield(options_generator,'rho')
        ff_name=[ff_name '_rho' num2str(options_generator.rho) '_lmq' num2str(options_generator.lambda_q)];         
    end
    ff_name=[ff_name '_algo_' pca_algorithm];
    
    save(fullfile(folder_exp,[ff_name '.mat']))
    disp(ff_name)
    disp(['reconstr_error:' num2str(mean(errors_reconstr))])

else
    errors_real=[];
    errors_batch_pca=[];
    errors_online=[];
    times_=[];
    ff_name='';
end

if ~compute_error_real
    errors_real=NaN;
end
if ~compute_error_batch
    errors_batch_pca=NaN;
end
if ~compute_error_online
    errors_online=NaN;
end

xlabel('Samples')
ylabel('Eigenspace Estimation Error Pq vs Pqhat')