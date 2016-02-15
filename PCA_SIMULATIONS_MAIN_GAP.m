clear all

options_simulations=struct;

options_simulations.outer_iter=1;
options_simulations.n0=0;
options_simulations.niter=40;
options_simulations.nstep_skip_EIGV_errors=128;

options_algorithm=struct;
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
options_generator.method='spiked_covariance_normalized';%'spiked_covariance';'brownian_motion';
%options_generator.rho=0.1;
options_generator.lambda_q=.95;
%%

%%
format compact
options_simulations.compute_error=1;

t = datetime('now');
folder_exp=['n0_' num2str(options_simulations.n0) '_niter' num2str(options_simulations.niter) '_' char(t)];
mkdir(folder_exp)
hold all
legends={};
options_generator.n=2048;
cm=hot(220);
counter=0;

for q=[2 4 16 64 256 512]
    for d=[4 16 64 256 1024]
        for rho=.1:.2:.9
            options_generator.q=q;
            options_generator.d=d;
            
            disp(rho)
            options_generator.rho=rho;
            options_algorithm=struct();
            options_algorithm.pca_algorithm='IPCA';
            options_algorithm.q=options_generator.q;
            options_algorithm.tol=1e-7;
            [errors_real,errors_batch_pca,errors_online,times_,fname]=Online_PCA_simulations(folder_exp,options_simulations,options_generator,options_algorithm);
            if ~isempty(errors_real)
                axs=[];
                counter=counter+1;
                cols=cm(counter,:);
                axs(1)=plot(median(errors_real,2),'d','Linewidth',2,'color',cols);
                axs(2)=plot(median(errors_batch_pca,2),'+','Linewidth',2,'color',cols);
                axs(3)=plot(median(errors_online,2),'*','Linewidth',2,'color',cols);
                legend(axs,{'real','batch pca','online'}, 'Interpreter', 'none')
                xlabel(fname)
                ylabel('Projection error')
                drawnow
            end
            options_algorithm=struct();
            options_algorithm.pca_algorithm='H_AH_NN_PCA';
            options_algorithm.q=options_generator.q;
            options_algorithm.update_method='ls';
            options_algorithm.tol=1e-5;
            [errors_real,errors_batch_pca,errors_online,times_,fname]=Online_PCA_simulations(folder_exp,options_simulations,options_generator,options_algorithm);
            if ~isempty(errors_real)
                axs=[];
                counter=counter+1;
                cols=cm(counter,:);
                axs(1)=plot(median(errors_real,2),'d','Linewidth',2,'color',cols);
                axs(2)=plot(median(errors_batch_pca,2),'+','Linewidth',2,'color',cols);
                axs(3)=plot(median(errors_online,2),'*','Linewidth',2,'color',cols);
                legend(axs,{'real','batch pca','online'}, 'Interpreter', 'none')
                xlabel(fname)
                ylabel('Projection error')
                drawnow
            end
        end
    end
end
%%
format compact
options_simulations.compute_error=0;
options_simulations.niter=10;

t = datetime('now');
folder_exp=['n0_' num2str(options_simulations.n0) '_niter' num2str(options_simulations.niter) '_' char(t)];
mkdir(folder_exp)
hold all
legends={};
options_generator.n=10000;
cm=hot(220);
counter=0;

for q=[128 256 512]
    for d=[65536]
        for rho=[.5]
            options_generator.q=q;
            options_generator.d=d;
            
            disp(rho)
            options_generator.rho=rho;
            options_algorithm=struct();
            options_algorithm.pca_algorithm='IPCA';
            options_algorithm.q=options_generator.q;
            options_algorithm.tol=1e-7;
            [errors_real,errors_batch_pca,errors_online,times_,fname]=Online_PCA_simulations(folder_exp,options_simulations,options_generator,options_algorithm);
            if ~isempty(errors_real)
                axs=[];
                counter=counter+1;
                cols=cm(counter,:);
                axs(1)=plot(median(errors_real,2),'d','Linewidth',2,'color',cols);
                axs(2)=plot(median(errors_batch_pca,2),'+','Linewidth',2,'color',cols);
                axs(3)=plot(median(errors_online,2),'*','Linewidth',2,'color',cols);
                legend(axs,{'real','batch pca','online'}, 'Interpreter', 'none')
                xlabel(fname)
                ylabel('Projection error')
                drawnow
            end
            options_algorithm=struct();
            options_algorithm.pca_algorithm='H_AH_NN_PCA';
            options_algorithm.q=options_generator.q;
            options_algorithm.update_method='ls';
            options_algorithm.tol=1e-5;
            [errors_real,errors_batch_pca,errors_online,times_,fname]=Online_PCA_simulations(folder_exp,options_simulations,options_generator,options_algorithm);
            if ~isempty(errors_real)
                axs=[];
                counter=counter+1;
                cols=cm(counter,:);
                axs(1)=plot(median(errors_real,2),'d','Linewidth',2,'color',cols);
                axs(2)=plot(median(errors_batch_pca,2),'+','Linewidth',2,'color',cols);
                axs(3)=plot(median(errors_online,2),'*','Linewidth',2,'color',cols);
                legend(axs,{'real','batch pca','online'}, 'Interpreter', 'none')
                xlabel(fname)
                ylabel('Projection error')
                drawnow
            end
        end
    end
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
%     [allerrors,alltimes,legends]=Online_PCA_simulations_1(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper);
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
% [allerrors,alltimes,legends]=Online_PCA_simulations_1(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,options_generator,errors_paper,errors_paper_half,times_paper,learning_rates);
% close all