function [M,W,Y,Ysq]=SEQ_ACT_SM_PCAFast(M,W,Y,Ysq,x,options)

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
    if ~isfield(options,'seq_act')
        seq_act=1;
    else
        seq_act=options.seq_act;
    end
    
    if ~isfield(options,'max_res')
        max_res=1e-1;
    else
        max_res=options.max_res;
    end

    if ~isfield(options,'tol')
        options.tol=1e-5;
    end

    if ~isfield(options,'update_method')
        options.update_method='ls';
    end
    
     if ~isfield(options,'lambda')
        options.lambda=0;
     end
    
    act=sum(Ysq>0);
    if act == options.q
        seq_act=0;
    end
    
    if seq_act        
        idx=1:act;
        Y=Y(idx);
        M=M(idx,idx);
        W=W(idx,:);
        Ysq=Ysq(idx);
        q=act;
    else     
        q=options.q;      
    end
    
    lambda=options.lambda;
    
    
    x=x';
    d=size(W,2);
    if isequal( options.update_method,'ls')
        % least square
%          disp(cond(eye(q)+M,'fro'))
%         pause
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
    
    if seq_act
        res=(sum(x.^2)-sum(Y.^2))/sum(x.^2);        
        if res>max_res
            act=act+1;
            Y(act,1)=sqrt(res*sum(x.^2));
            Ysq(act,1)=0;
            M(act,act)=0;
            W(act,:)=0;            
            %disp('** adding component **')
        end
    end
    
    Ysq = Ysq + Y.^2;
    
    % Update weights
    Y_tmp=(Y./Ysq);
    Y_tmp_sq=Y.^2./Ysq;

    % Update weights
    %W = W + Y_tmp*x' - W.*repmat(Y_tmp_sq,[1 d]);
     W = W +  bsxfun(@times,Y_tmp,x') - bsxfun(@times,W,Y_tmp_sq);


    if isnan(sum(W(:)))
        W(isnan(W)) = 0;
    end

    M = M + (1+lambda)*bsxfun(@times,Y_tmp,Y') - bsxfun(@times,M,Y_tmp_sq);
   % M = M + Y_tmp*Y' - M.*repmat(Y_tmp_sq,[1 q]);       
    M(isnan(M)) = 0;
    %stupid comment
    M(1:q+1:end)=0;
    
    if seq_act        
        M(act+1:options.q,act:options.q)=0;
        Y(act+1:options.q,1)=0;
        Ysq(act+1:options.q,1)=0;
        W(act+1:options.q,:)=0;
    end
end
