function [M,W,Y]=H_AH_NN_PCAFast(M,W,Y,x,options)
    % function implementing Hebbian Antihebbian learning for subspace learning
    % [M,W,Y]=H_AH_NN_PCAFast(M,W,Y,x,options)
    % W: forward connection matrix
    % M: lateral connection matrix
    % Y: 
    % options.update_method: method to perform the update of Y variable. 'ls' (least-square),
    % 'coord_desc' (coordinate descent), or 'mat_mult updated' (update all coordinated simultaneously)
    % options.tol: tolerance on convergence.
    % reference: Pehlevan et al, Neural Computation, 2015   
    % gamma=1./Ysq;
    
    
    if ~isfield(options,'tol')
        options.tol=1e-5;
    end

    if ~isfield(options,'update_method')
        options.update_method='ls';
    end
    
    q=options.q;      
    gamma=options.gamma;    
    
    
    x=x';
    d=size(W,2);
    if isequal( options.update_method,'ls')
        % least square
        Y=(eye(q)+M)\(W*x);
    else
        er = Inf;
        iter_num=0;
%         Y=rand(size(W,1),1);
        while er > options.tol && iter_num <= options.mat_iter
            iter_num=iter_num+1;
            Yprev = Y;
            switch options.update_method
                case 'coord_desc'
                    % Coordinate descent step until convergence
                    for i = 1:q
                        Y(i) = W(i,:)*x - M(i,:)*Y;
                    end
                case 'mat_mult'
                    % matrix multiplication until convergence
                    Y= W*x - M*Y;
                otherwise
                    error('Unrecognized method')

            end
            er = max(abs(Y-Yprev)./abs(Yprev));
        end
    end
    
    % Update weights
    Y_tmp=(Y.*gamma);
    Y_tmp_sq=Y.^2.*gamma;

    % Update weights
    %W = W + Y_tmp*x' - W.*repmat(Y_tmp_sq,[1 d]);
    W = W +  bsxfun(@times,Y_tmp,x') - bsxfun(@times,W,Y_tmp_sq);

    if isnan(sum(W(:)))
        W(isnan(W)) = 0;
    end

    M = M + bsxfun(@times,Y_tmp,Y') - M.*repmat(Y_tmp_sq,[1 q]);
   % M = M + Y_tmp*Y' - bsxfun(@times,M,Y_tmp_sq);
    M(isnan(M)) = 0;
    %stupid comment
    M(1:q+1:end)=0;
    
end
