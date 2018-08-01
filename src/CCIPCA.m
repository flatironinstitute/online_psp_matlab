classdef CCIPCA < handle
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
      ell
      v
      
   end
   
   methods
      function obj = CCIPCA(k, d, Uhat0, lambda0, ell0)
          obj.k = k;
          obj.d = d;
          if isempty(Uhat0)   
              obj.Uhat = randn(k,d)/d;
          else
              obj.Uhat = Uhat0;
          end
          if isempty(lambda0)
              obj.lambda_ = abs(randn(k,1)/sqrt(k));
          else
              obj.lambda_ = lambda0;
          end
          if isempty(ell0)
              obj.ell = 2;
          else
              obj.ell = ell0;
          end
          obj.t = 1;
          obj.v = zeros(d,1);
      end
      
      function fit_next(obj,x)
        old_wt = max(1, obj.t - obj.ell) / (obj.t + 1);
        for i = 1:obj.k
            obj.v = old_wt * obj.lambda_(i) * obj.Uhat(i, :) + (1 - old_wt) * x * obj.Uhat(i, :)' .* x;
            obj.lambda_(i) = norm(obj.v);
            obj.Uhat(i, :) = obj.v / obj.lambda_(i);
            x = x - x*obj.Uhat(i, :)'.*obj.Uhat(i, :);            
            obj.t = obj.t + 1;
        end
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
              
          components = obj.Uhat';
          
          if orthogonalize
                [components, ~ ] = qr(components,0);
          end  
          
      end
      
   end
end
