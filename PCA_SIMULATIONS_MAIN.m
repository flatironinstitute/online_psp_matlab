clear all
profile off
profile on
niter=50; % was 40
n0 = 250; % number of sample paths used for initialization
qreal = 5;% was 5	 % number of PCs to çcompute
q=1*qreal;
q_each=[5];
num_samples_each=[1000];
d_each=[10,100,1000];

if n0==50
    learning_rates=[10,1,.1]/5; 
else 
    learning_rates=[10,1,.1];
end

nstep_skip_EIGV_errors=50;
q_mult=1;
close all
%%
options=struct();
test_method='H_AH_NN_PCA';
options.tol=1e-5;
%options.update_method='coord_desc';
errors_paper=[.011 .007 .007];
errors_paper_half=[0.02 .015 .015];
times_paper=[.17  .11  .17 ];
[allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,q_mult,n0,niter,nstep_skip_EIGV_errors,test_method,options,errors_paper,errors_paper_half,times_paper);
close all
%%
options=struct();
test_method='IPCA';
options.tol=1e-7;
errors_paper=[.011 .007 .007];
errors_paper_half=[0.02 .015 .015];

times_paper=[.17  .11  .17 ];
[allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,q_mult,n0,niter,nstep_skip_EIGV_errors,test_method,options,errors_paper,errors_paper_half,times_paper);
close all
%% !!!! LEARNING RATE MANAULLY SET IN THE FILE!!!
options=struct();
test_method='SGA'; 
options.method=test_method;
options.do_sort=1;
errors_paper=[ .025 .014 .016];
errors_paper_half=[0.031 .02 .021];

times_paper=[0.09 .16 .12];
[allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,q_mult,n0,niter,nstep_skip_EIGV_errors,test_method,options,errors_paper,errors_paper_half,times_paper,learning_rates);
close all
%% !!!! LEARNING RATE MANAULLY SET IN THE FILE

options=struct();
test_method='GHA'; 
options.method=test_method;
options.do_sort=1;
errors_paper=[ .024 .014 .016];
errors_paper_half=[0.03 .02 .023 ];

times_paper=[0.06 .08 .09 ];
[allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,q_mult,n0,niter,nstep_skip_EIGV_errors,test_method,options,errors_paper,errors_paper_half,times_paper,learning_rates);
close all