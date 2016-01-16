function [allerrors,alltimes,legends]=Online_PCA_simulations(q_each,num_samples_each,d_each,n0,niter,nstep_skip_EIGV_errors,test_method,options,errors_paper,errors_paper_half,times_paper,learning_rates)
%%
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
                [I,J] = ind2sub([d,d],1:d^2);
                C=zeros(d);
                C(:)=min(I,J)/d;
                [Q,L,~]=eigs(C,q);
                Pq=Q(:,1:q)*Q(:,1:q)';
                % n = 500; % number of sample paths
                % d = 10;	 % number of observation points
                errors=zeros(n,niter);
                times_=zeros(n,niter);
                for ll=1:niter
                    disp(ll)
                    x1 = normrnd(0,1/sqrt(d),n,d);
                    x= cumsum(x1,2);
                    
                    %             [U1,S1,V1]=svds(x(1:n0,:),q);
                    %             eigval1=diag(S1*S1)/(n0-1);
                    %             subplot(3,1,1)
                    %             imagesc(V1)
                    %             subplot(3,1,2)
                    %             imagesc(COEFF)
                    %             subplot(3,1,3)
                    %             imagesc(Q)
                    
                    [eigvect,~,eigval]=pca(x(1:(n0+1),:),'NumComponents',q);
                    %             [egve,egvl]=svd(x(1:n0,:));
                    %             [V,D,W] = svd(x(1:n0,:)'*x(1:n0,:)/(n0-1));
                    %             %[COEFF,SCORE,latent]=pca(x(1:n0,:));
                    %             eigval=diag(D);
                    %             values=diag(D);
                    %             eigvect=V;
                    
                    
                    values=eigval(1:q);
                    vectors=eigvect(:,1:q);
                    
                    if isequal('H_AH_NN_PCA',test_method)
                        M = 0.0*rand(q,q);
                        M(1:q+1:end)=0; % set diagonals to zero
                        W = vectors';
                        Y = zeros(q,1);
                        %Ysq=max(10*ones(size(W,1),1),sum((W*x(1:n0,:)').^2,2));
                        Ysq=10*ones(size(W,1),1)+sum((W*x(1:n0,:)').^2,2);
                    end
                    
                    
                    Pq_hat=vectors(:,1:q)*vectors(:,1:q)';
                    errors(n0,ll)=(norm(Pq_hat-Pq,'fro')^2)/(norm(Pq,'fro')^2);
                    
                    
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
                            Pq_hat=vectors(:,1:q)*vectors(:,1:q)';
                            errors(i,ll)=(norm(Pq_hat-Pq,'fro')^2)/(norm(Pq,'fro')^2);
                        end
                        %         errors(i,ll)=2*(1-trace(Pq_hat*Pq)/q);
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
