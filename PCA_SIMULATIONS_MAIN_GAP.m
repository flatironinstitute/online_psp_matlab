clear all
profile off
profile on
niter=40; % was 40
n0 = 32; % number of sample paths used for initialization
q_each=16;
num_samples_each=4096;
d_each=[256];

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
%%
options_generator=struct;
options_generator.method='spiked_covariance_normalized';
% options_generator.rho=0.1;
options_generator.lambda_q=.95;
%%
errors_paper=[];
errors_paper_half=[];
times_paper=[];
for rho=.8%0:.1:1
    disp(rho)
    options_generator.rho=rho;    
    options=struct();
    test_method='IPCA';
    options.tol=1e-7;        
    [~,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
    options=struct();
    test_method='H_AH_NN_PCA';
    options.update_method='ls';
    options.tol=1e-5;
    [allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
%     uiwait
%     cla
end
%%
errors_paper=[];
errors_paper_half=[];
times_paper=[];
for rho=0:.1:1
    disp(rho)
    options_generator.rho=rho;    
    options=struct();
    test_method='SGA';
    options.method=test_method;
    options.tol=1e-7;  
    options.do_sort=1;
    [allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
    %     uiwait
%     cla
end


%% !!!! LEARNING RATE MANUALLY SET IN THE FILE!!!
% for rho=0:.1:1
%     disp(rho)
%     options_generator.rho=rho;    
%     
%     options=struct();
%     test_method='SGA';
%     options.method=test_method;
%     options.do_sort=1;
%     errors_paper=[ .025 .014 .016];
%     errors_paper_half=[0.031 .02 .021];
%     
%     times_paper=[0.09 .16 .12];
%     [allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
% end
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