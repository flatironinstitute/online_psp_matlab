%%Use ONLINE PCA with Hebbian Anti-Hebbian Learning rule (CEnzig et. al,2015 Neural Computation)
clear all
load fisheriris
x = meas';

n_init_PCA = 0; % number of samples used for initialization
d=size(x,1); % input dimensionality
k=2; % output dimensionality
n=size(x,2);
scrambled = randperm(n);
x = x(:,scrambled);
%% BATCH PCA
[coeff,score,pcvar] = pca(x','NumComponents',k);
%% ONLINE PCA
[x,eig_vect,eig_val] = standardize_data(x,coeff,pcvar);
Uhat0 = bsxfun(@times,x( :,1:k), 1./sqrt(sum(x(:,1:k).^2,1)))';
scal = 100;
Minv0 = eye(k) * scal;
Uhat0 = Uhat0 / scal;
fsm = FSM(k, d, [], Minv0, Uhat0, []);
Y = zeros(k,n);
for i = 1:n
    if mod(i,50) == 1
        disp(i)
    end
    Y(:,i) = fsm.fit_next(x(:,i)');
end
%% PLOT RESULT BATCH
figure
scatterhist(score(:,1),score(:,2), 'group',species(scrambled))
%% PLOT RESULTS ONLINE
figure
scatterhist(Y(1,:),Y(2,:), 'group',species(scrambled(n_init_PCA+1:end)))
