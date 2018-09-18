classdef IPCA < handle
    %{
    """
    Candid covariance Incremental PCA
    Parameters:
    ====================
    X             -- Numpy array of size d-by-n, where each column corresponds to one observation
    q             -- Dimension of PCA subspace to learn, must satisfy 1 <= q <= d
    Uhat0         -- Initial guess for the eigenspace matrix U, must be of size d-by-q
    lambda0       -- Initial guess for the eigenvalues vector lambda_, must be of size q
    ell           -- Amnesiac parameter (see reference)

    Methods:
    ====================
    fit()
    %}
   properties
      Uhat
      t
      lambda_
      k
      d
      f
      tol
      
   end
   
   methods
      function obj = IPCA(k, d, Uhat0, lambda0, f0, tol0)
          obj.k = k;
          obj.d = d;
          if isempty(Uhat0)   
              obj.Uhat = randn(d,k)/d;
          else
              obj.Uhat = Uhat0';
          end
          if isempty(lambda0)
              obj.lambda_ = abs(randn(k,1)/sqrt(k));
          else
              obj.lambda_ = lambda0;
          end
          obj.t = 1;
          if isempty(f0)
              obj.f = 1.0 / obj.t;
          else
              obj.f = f0;
          end
          if isempty(tol0)
              obj.tol = 1e-7;
          else
              obj.tol = tol0;
          end
          
      end
      
      function fit_next(obj,x)
          x = x';
          obj.t = obj.t + 1;
          obj.f = 1.0 / obj.t;          
          obj.lambda_ = (1 -  obj.f) *  obj.lambda_;
          x = sqrt(obj.f) * x;
          % Project X into current estimate and check residual error
          Uhatx =  obj.Uhat' * x;
          x = x -  obj.Uhat * Uhatx;
          normx = sqrt(x'*x);  % norm(x)
          
          % TODO: fix this atleast_2d for efficiency
          if (normx >= obj.tol)
              obj.lambda_ = cat(1,obj.lambda_,0);
              Uhatx = cat(1, Uhatx, normx);
              obj.Uhat = cat(2,obj.Uhat, x / normx);
          end
          M = diag( obj.lambda_) + Uhatx * Uhatx';
          [V, s] = eig(M);
          s = diag(s);
          if (normx >= obj.tol)
            obj.lambda_ = s(end:-1:2);
            V = V(:, end:-1:2);       
          end
          obj.Uhat =  obj.Uhat*V;
       
      end
      
      function components = get_components(obj, orthogonalize)       
        %{       
        Extract components from object
         Parameters
         ---------
         orthogonalize: bool
             whether to orthogonalize when computing the error
 
         Returns
         -------
         components: ndarray 
         %}
          
          if isempty(orthogonalize)
              orthogonalize = 1;
          end
              
          components = obj.Uhat;
          
          if orthogonalize
                [components, ~ ] = qr(components,0);
          end  
          
      end
      
   end
end


