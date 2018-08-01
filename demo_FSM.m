%%Use ONLINE PCA with Hebbian Anti-Hebbian Learning rule (Pehlevan et. al,2015 Neural Computation)
clear all
d=256; % input dimensionality
k=16; % output dimensionality
n=6000;
method = 'spiked_covariance_normalized';
options_generator=struct;
options_generator.rho=0.01;
options_generator.lambda_q=1;
[x,eig_vect,eig_val] = low_rank_rnd_vector(d,k,n,method,options_generator);
[x,eig_vect,eig_val] = standardize_data(x,eig_vect,eig_val);
disp('generated')
%%
errors=zeros(n,1)*NaN;
Uhat0 = bsxfun(@times,x( :,1:k), 1./sqrt(sum(x(:,1:k).^2,1)))';
scal = 100;
Minv0 = eye(k) * scal;
Uhat0 = Uhat0 / scal;
fsm = FSM(k, d, [], Minv0, Uhat0, []);

for i = 1:n
    if mod(i,50) == 1
        errors(i,:) = compute_projection_error(eig_vect, fsm.get_components([]));
        disp(i)
    end
    fsm.fit_next(x(:,i)');
end
loglog(1:n,errors,'.')