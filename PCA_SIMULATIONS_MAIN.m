clear all
profile off
profile on
niter=10; % was 40
n0 = 0; % number of sample paths used for initialization
q_each=[2 4 16 64 256];% 512];
num_samples_each=2048;%8192;
d_each=[4 16 64 256 1024]% 4096];

nstep_skip_EIGV_errors=50;
close all
%%
% options_generator=struct;
% options_generator.method='brownian_motion';
%%
% options_generator=struct;
% options_generator.method='spiked_covariance';
% options_generator.rho=0.1;
% options_generator.gap=.5;
% options_generator.slope=.05;
% options_generator.compute_eig=1;
%%
options_generator=struct;
options_generator.method='spiked_covariance_normalized';
options_generator.rho=0.8;
options_generator.lambda_q=.95;
%%
options=struct();
test_method='H_AH_NN_PCA';
options.tol=1e-5;
% options.update_method='mat_mult';
% options.mat_iter=10;
options.update_method='ls';
% errors_paper=[.011 .007 .007];
% errors_paper_half=[0.02 .015 .015];
errors_paper=[];
errors_paper_half=[];
times_paper=[.17  .11  .17 ];
[allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
close all
%%
options=struct();
test_method='IPCA';
options.tol=1e-7;
errors_paper=[.011 .007 .007];
errors_paper_half=[0.02 .015 .015];

times_paper=[.17  .11  .17 ];
[allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
close all
%% !!!! LEARNING RATE MANAULLY SET IN THE FILE!!!
% options=struct();
% test_method='SGA'; 
% options.method=test_method;
% options.do_sort=1;
% errors_paper=[ .025 .014 .016];
% errors_paper_half=[0.031 .02 .021];
% 
% times_paper=[0.09 .16 .12];
% [allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
% close all
% %% !!!! LEARNING RATE MANAULLY SET IN THE FILE
% 
% options=struct();
% test_method='GHA'; 
% options.method=test_method;
% options.do_sort=1;
% errors_paper=[ .024 .014 .016];
% errors_paper_half=[0.03 .02 .023 ];
% 
% times_paper=[0.06 .08 .09 ];
% [allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper,learning_rates);
% close all