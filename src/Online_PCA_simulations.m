function [allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,errors_paper,errors_paper_half,times_paper,learning_rates)
%%
method_random='brownian_motion';
options=struct;

close all
set(0, 'DefaulttextInterpreter', 'none')
subplot(2,2,1)
t = datetime('now');
folder_exp=[test_method '_n0' num2str(n0) '_' char(t)];
mkdir(folder_exp)
allerrors=[];
allerrors_half=[];
alltimes={};
legends={};
for q=q_each
    options.q=q;
    for n=num_samples_each
        counter=0;
        for d=d_each
            counter=counter+1;            
            if d>=q     
                if isequal(test_method,'SGA') || isequal(test_method,'GHA')
                    learning_rate=learning_rates(counter);
                end
                hold all;
                    
                errors=zeros(n,niter);
                times_=zeros(n,niter);
                for ll=1:niter
                    disp(ll)
                    %generate random samples
                    if isequal(method_random,'brownian_motion') && ll>1 % eigenvalues are only computed once since the covariance matrix is always the same
                        options.compute_eig=0;
                        [x,~,~] = low_rank_rnd_vector(d,q,n,method_random,options);
                    else
                        options.compute_eig=1;
                        [x,eig_vect,eig_val] = low_rank_rnd_vector(d,q,n,method_random,options);                        
                    end
                    
                    x1 = normrnd(0,1/sqrt(d),n,d);
                    x= cumsum(x1,2);
                    
                    [eigvect_init,~,eigval_init]=pca(x(1:(n0+1),:),'NumComponents',q);

                    values=eigval_init;
                    vectors=eigvect_init;
                    
                    if isequal('H_AH_NN_PCA',test_method)
                        M = zeros(q,q);
%                         M(1:q+1:end)=0; % set diagonals to zero
                        W = vectors';
                        Y = zeros(q,1);
                        %Ysq=max(10*ones(size(W,1),1),sum((W*x(1:n0,:)').^2,2));
                        Ysq=0*10*ones(size(W,1),1)+sum((W*x(1:n0,:)').^2,2);
                    end
                    
                    
                   
                    errors(n0,ll)=compute_reconstruction_error(eig_vect,vectors);
                    
                    
                    for i = (n0+1):n
                        tic
                        %disp(i)
                        switch test_method
                            case 'SGA'
                                options.gamma=learning_rate/i;
                                [values, vectors] = sgaPCAFast(values, vectors, x(i,:),options);
                            case 'GHA'
                                options.gamma=learning_rate/i;
                                [values, vectors] = sgaPCAFast(values, vectors, x(i,:),options);
                            case 'IPCA'
                                options.n=i-1;
                                [values, vectors] = incrPCA_fast(values, vectors,  x(i,:), options);
                            case 'H_AH_NN_PCA'
                                options.gamma=1./Ysq;
                                [M,W,Y] = H_AH_NN_PCAFast(M,W,Y, x(i,:),options);
                                Ysq = Ysq + Y.^2;
                                
                        end
                        times_(i,ll)=toc;
                        if mod(i,nstep_skip_EIGV_errors) ==  0 || i==n || i == round(n/2) 
                            if isequal('H_AH_NN_PCA',test_method)
                                %disp(['iteration:' num2str(i) ', computing vectors...'])
                                vectors = (pinv(diag(ones(q,1))+M(1:q,1:q))*W(1:q,:))';
                            end                         
                            errors(i,ll)=compute_reconstruction_error(eig_vect,vectors);
                        end
                    end
                    %             plot(errors(:,ll))
                    %             drawnow
                end
                
                allerrors=[allerrors; errors(end,:)];
                allerrors_half=[allerrors_half; errors(round(n/2),:)]; %remember that most are zeros!!
                alltimes=[alltimes {times_(times_~=0)}];
                plot(median(errors,2),'d','Linewidth',2)
                drawnow
                save(fullfile(folder_exp,['n' num2str(n) '_d' num2str(d) '_q' num2str(q)]))
                legends=[legends ['n' num2str(n) '_d' num2str(d) '_q' num2str(q)]];
                
            end
        end
    end
end
legend(legends, 'Interpreter', 'none')
xlabel('Samples')
ylabel('Eigenspace Estimation Error Pq vs Pqhat')
%%
subplot(2,2,2)
errorbar(cellfun(@median,alltimes)*1000,cellfun(@iqr,alltimes)*1000)
hold on
plot(times_paper,'r*')
ax=gca;
ax.XTick=[1:numel(legends)];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
ylabel('Time per iteration (ms)')

subplot(2,2,3)
hold on
errorbar(median(allerrors_half'),iqr(allerrors_half'),'c')
% boxplot(allerrors_half','whisker',100)
plot(errors_paper_half,'k*')
ax=gca;
ax.XTick=[1:numel(legends)];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
ylabel('Half Samples Eigenspace Estimation Error Pq vs Pqhat')

subplot(2,2,4)
hold on
errorbar(median(allerrors'),iqr(allerrors'))
% boxplot(allerrors','whisker',100)

plot(errors_paper,'k*')
ax=gca;
ax.XTick=[1:numel(legends)];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
ylabel('Eigenspace Estimation Error Pq vs Pqhat')

savefig(fullfile(folder_exp,'times_errors.fig'))
