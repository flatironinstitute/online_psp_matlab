function [M,W,Y]=minMaxPCA(M,W,Y,x,eta,tau,options)

    % function implementing Hebbian Anti-hebbian networks for subspace learning
    % [M,W,Y]=H_AH_NN_PCAFast(M,W,Y,x,options)
    % W: forward connection matrix
    % M: lateral connection matrix
    % Y: output space
    % options.update_method: method to perform the update of Y variable. 'ls' (least-square),
    % 'coord_desc' (coordinate descent), or 'mat_mult updated' (update all coordinated simultaneously)
    % options.tol: tolerance on convergence.
    % options.lambda: factor influencing decorrelation (0 no decorrelatoin applied, see NIPS 2015)
    % reference: Pehlevan et al, Neural Computation, 2015. Pehlevan et al, NIPS, 2015   
    % gamma=1./Ysq;
    
    
    if ~isfield(options,'tol')
        options.tol=1e-5;
    end

    if ~isfield(options,'update_method')
        options.update_method='ls';
    end
    
     if ~isfield(options,'lambda')
        options.lambda=0;
     end
    
    
    q=options.q;      
    lambda=options.lambda;
    
    x=x';
    d=size(W,2);
    if isequal( options.update_method,'ls')
        % least square
%          disp(cond(eye(q)+M,'fro'))
%         pause
        Y=(M)\(W*x);
    else
        er = Inf;
        iter_num=0;
%         Y=rand(size(W,1),1);
        while er > options.tol && iter_num <= options.mat_iter
            iter_num=iter_num+1;
            error('Need to decide how to initialize Y!')
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
    
%     Ysq = Ysq + Y.^2;
%     
%     % Update weights
%     Y_tmp=(Y./Ysq);
%     Y_tmp_sq=Y.^2./Ysq;

    % Update weights
    %W = W + Y_tmp*x' - W.*repmat(Y_tmp_sq,[1 d]);
     W = W + 2*eta*(Y*x'-W);
     M = M + eta/tau*(Y*Y'-M);

    if isnan(sum(W(:)))
        W(isnan(W)) = 0;
    end
    
end
