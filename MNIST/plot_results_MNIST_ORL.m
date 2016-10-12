%% display results simulation
clear all
files=dir('*.mat')
figure
hold all
legends={};
cm=colormap(hot(numel(files)+10));
d_s=[];
q_s=[];
err_s=[];
errhalf_s=[];
for f=1:numel(files)
   disp(files(f).name)
   load(files(f).name,'d','q','errors_real','errors_batch_pca')
   legends{f}=files(f).name;
   errors=errors_batch_pca;
   plot((median(errors,2)),'d','Linewidth',2,'Color',cm(f+5,:))
   err_s=[ err_s errors(end,:)'];
   d_s=[d_s d];
   q_s=[q_s q];
   legend(legends{:},'Interpreter', 'none')   
   xlabel(files(f).name)
   ylim([0 2])
%    pause
%    cla
%     input('')

end
legend(legends, 'Interpreter', 'none')
xlabel('Samples')
ylabel('Eigenspace Estimation Error pca vs batch')
%%
clear all
files=dir('*.mat');
% dirFlags = [files.isdir];
% files=files(dirFlags);
files = struct2cell(files);
files = files(1,:);
%%
files=uipickfiles;
%% ********************************** PLOT FOR GAPS *****************************************
clear all
files_to_analize = uipickfiles()
%%
files_to_analize=dir('*')
files_to_analize(1:2)=[];
dirFlags = [files_to_analize.isdir];
files_to_analize=files_to_analize(dirFlags);
files_to_analize = struct2cell(files_to_analize);
files_to_analize = files_to_analize(1,:);
%%
files_to_analize=dir('*.mat')
% dirFlags = [files.isdir];
% files=files(dirFlags);
files_to_analize = struct2cell(files_to_analize);
files_to_analize = files_to_analize(1,:);
%% PLOT TOP EIGENVALUES
files_to_analize = uipickfiles();
% d1=112;
% d2=92;
d1=28;
d2=28;
% d1=32;
% d2=32;
num_comps=5;
for mm=1
    load(files_to_analize{1})
    F=(pinv(diag(ones(q,1))+M(1:q,1:q))*W(1:q,:))';
    F = orth(F);
    pwr=sum((F'*x').^2,2);
    [vls,idx]=sort(pwr,'descend');
    plot(pwr(idx))

    for cc=1:num_comps
       subplot(3,5,cc)
       imagesc(reshape(eig_vect_batch_pca(:,cc),[d1,d2]))
       axis image
       colormap gray
        axis off

       title(['PCA comp:' num2str(cc)])
    end
    
    for cc=1:num_comps
       subplot(3,5,cc+2*num_comps)
       imagesc(reshape(F(:,idx(cc)),[d1,d2]))
       axis image
       axis off
       colormap gray
       title(['OSM comp:' num2str(cc)])
    end
    load(files_to_analize{2})
     
    for cc=1:num_comps
       subplot(3,5,cc+1*num_comps)
       imagesc(reshape(vectors(:,cc),[d1,d2]))
       axis image
       colormap gray
              axis off

       title(['IPCA comp:' num2str(cc)])
    end
end
saveas(gcf,'first5comp.fig')
% exportfig(gcf,'YALEFirst5.eps','color','cmyk','fontsize',1,'width',10,'height',10)
%% PLOT ERROR IN FUNCTION OF TRIALS
files_to_analize = uipickfiles();
hold all
load(files_to_analize{1})
err2=nanmedian(errors_batch_pca');
load(files_to_analize{2})
err1=nanmedian(errors_batch_pca');
idxs=find(~isnan(err2-err1));
plot(idxs,err2(idxs)-err1(idxs),'-d')

ylabel('\Delta Projection error')
xlabel('Samples')
legend('MNIST','YALE','ATT','ORL','SPIKED-COV')
%%
files_to_analize = uipickfiles();
hold all
load(files_to_analize{1})
err2=nanmedian(errors_batch_pca');
load(files_to_analize{2})
err1=nanmedian(errors_batch_pca');
idxs=find(~isnan(err2-err1));


plot(idxs,err2(idxs),'--*k')
plot(idxs,err1(idxs),'--*r')

ylabel('Projection error')
xlabel('Samples')
legend('OSM YALE','IPCA YALE','OSM ATT','IPCA ATT')
%% PLOT SUMMARY OF FILES
clear all
files=dir('*.mat');
files = struct2cell(files);
files = files(1,:);
load_times=0;

times_tot=[];
err_real=[];
err_batch=[];
err_online=[];
err_reconstr=[];
d_s=[];
q_s=[];
rho_s=[];
methods_={};
for f=1:numel(files)
   disp(f)  
   load(files{f},'options_algorithm','options_generator','options_simulations','d','q','times_','errors_real','errors_batch_pca','errors_online','errors_reconstr')
   if ~isfield(options_generator,'rho')
        warning('setting rho to 0 because not existing field!')
        options_generator.rho=0;
   end
   test_method=options_algorithm.pca_algorithm;  
   if load_times
       times_=diff(times_);
       times_iter=times_(:)';
    %    times_iter=nanmean(times_,1); 
       numIter=numel(times_iter);
       times_tot=[times_tot times_iter];
   else
       numIter=size(errors_batch_pca,2);
   end

   idx_not_nan=find(~isnan(nanmedian(errors_batch_pca,2)));
   if ~load_times && ~isempty(idx_not_nan)       
       err_real=[err_real errors_real(idx_not_nan(end),:)];
       err_batch=[err_batch errors_batch_pca(idx_not_nan(end),:)];
       err_online=[err_online nanmean(errors_online(idx_not_nan,:),1)];  
       err_reconstr=[err_reconstr errors_reconstr(idx_not_nan(end),:)];
   end
   d_s=[d_s repmat(d,1,numIter)];
   q_s=[q_s repmat(q,1,numIter)];
   rho_s=[rho_s repmat(options_generator.rho,1,numIter)];
   newm={};
   for ll=1:numIter
       newm{ll}=test_method;
   end
   methods_=[methods_ newm];
end


%% error plot
figure('name','5_percentile')

is_proj_error=0;

if is_proj_error
    error=err_batch;
else
    error=err_reconstr;
end

jj=0;



cm1=hot(8);
cm2=(gray(10));
cm3=(autumn(10));
cm4=(lines(10));

stats_={'nanmean','sem'};
% stats_={@(X) quantile(X,.5),@(X) nan};

col_var=d_s;
col_var2=rho_s;
for cv=unique(col_var)%[400 784 1024]
    for cv2=unique(col_var2)
        legend off
        if ~isempty(find(col_var==cv & col_var2==cv2,1))
            jj=jj+1;
            subplot(1,1,jj)
            
            idx=find(col_var==cv & col_var2==cv2 & strcmp(methods_,'H_AH_NN_PCA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_h,ma_h]=grpstats(error(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.01,size(xax)), me_h,ma_h,'o-','MarkerSize',7,'MarkerFaceColor',cm1(5,:),'color',cm1(5,:));
            
            
            hold on
            idx=find(col_var==cv & col_var2==cv2 & strcmp(methods_,'IPCA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_i,ma_i]=grpstats(error(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.001,size(xax)), me_i,ma_i,'o-','MarkerSize',7,'MarkerFaceColor',cm2(5,:),'color',cm2(5,:));
            
             idx=find(col_var==cv & col_var2==cv2 & strcmp(methods_,'SEQ_SIM_PCA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_i,ma_i]=grpstats(error(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.001,size(xax)), me_i,ma_i,'o-','MarkerSize',7,'MarkerFaceColor',cm3(10,:),'color',cm3(10,:));
             
            idx=find(col_var==cv & col_var2==cv2 & strcmp(methods_,'SGA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_i,ma_i]=grpstats(error(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.001,size(xax)), me_i,ma_i,'o-','MarkerSize',7,'MarkerFaceColor',cm3(5,:),'color',cm3(5,:));
             
            idx=find(col_var==cv & col_var2==cv2 & strcmp(methods_,'CCIPCA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_i,ma_i]=grpstats(error(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.001,size(xax)), me_i,ma_i,'o-','MarkerSize',7,'MarkerFaceColor',cm4(5,:),'color',cm4(5,:));
            
            
%             xlim([3 300])
%             ylim([.001 1])
            axis tight
%              set(gca,'xscale','log')
            set(gca,'yscale','log')
            %axis tight
            xlabel('q')
            if is_proj_error
            ylabel('Projection Error')
            else
             ylabel('Reconstruction Error')   
            end
            %title(['d=' num2str(cv) ' q=' num2str(cv2)])
%             L = get(gca,'XLim');
%             set(gca,'XTick',xax)

            legend('OSM','IPCA')
            box off
        end
    end
end
set(gcf,'defaulttextinterpreter','none')
if is_proj_error
    saveas(gcf,'projErr.fig')
    exportfig(gcf,'projErr.eps','color','cmyk','fontsize',2,'width',7,'height',12)
else
    saveas(gcf,'reconstrErr.fig')
    exportfig(gcf,'reconstrErr.eps','color','cmyk','fontsize',2,'width',7,'height',12)

end
%% time plot
figure
jj=0;
ttime=times_tot;

cm1=hot(8);
cm2=(gray(10));
cm3=(autumn(10));
% stats_={'nanmean',@(X) mad(X,1),@(X) quantile(X,.25),@(X) quantile(X,.75)};
stats_={'nanmean',@(X) iqr(X),@(X) quantile(X,.25),@(X) quantile(X,.75)};

col_var=d_s;
for cv=unique(col_var)
    if ~isempty(find(col_var==cv,1))
        if ~isempty(find(col_var==cv ,1))
            jj=jj+1;
            subplot(1,1,jj)
            
            idx=find(col_var==cv & strcmp(methods_,'H_AH_NN_PCA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_h,ma_h,ql_h,qh_h]=grpstats(ttime(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.01,size(xax)), me_h,ma_h,'o-','MarkerSize',6,'MarkerFaceColor',cm1(5,:),'color',cm1(5,:));
            
            hold on
            idx=find(col_var==cv & strcmp(methods_,'IPCA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_i,ma_i,ql_i,qh_i]=grpstats(ttime(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.01,size(xax)), me_i,ma_i,'o-','MarkerSize',6,'MarkerFaceColor',cm2(5,:),'color',cm2(5,:));
            
            idx=find(col_var==cv & strcmp(methods_,'SGA'));
            xvar=q_s(idx);
            xax=unique(xvar);
            [me_i,ma_i,ql_i,qh_i]=grpstats(ttime(idx),xvar,stats_);
            errorbar(xax+normrnd(0,.01,size(xax)), me_i,ma_i,'o-','MarkerSize',6,'MarkerFaceColor',cm3(5,:),'color',cm3(5,:));
            
            
            
            legend('H_AH_NN_PCA','IPCA','SGA')
            set(gca,'yscale','log')
            set(gca,'xscale','log')
            xlabel('q')
            ylabel('time (ms)')
            L = get(gca,'XLim');
            set(gca,'XTick',xax)
            %ylim([2*1e-4 2])
            axis square
            title(['d=' num2str(cv)])
        end
    end
end
%% error plot online
clear all
files_to_analize = uipickfiles()
%%
figure('name','error,_online')

for jj=1:numel(files_to_analize)
    load(files_to_analize{jj})
    idx=find(~isnan(quantile(errors_real',.05)));
    if ~isfield(options_generator,'rho')
        warning('setting rho to 0 because not existing field!')
        options_generator.rho=0;
   end
    subplot(4,3,jj)
    
    hold all
    plot(idx,quantile(errors_real(idx,:)',.05),'-o')
    plot(idx,quantile(errors_online(idx,:)',.05),'-o')
    plot(idx,quantile(errors_batch_pca(idx,:)',.05),'-o')
    xlabel('Samples')
    ylabel('Projection error')
    title(['d=' num2str(d) ',q=' num2str(q) ',algo=' pca_algorithm ',rho=' num2str(options_generator.rho)])
    axis tight
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca, 'box', 'off')
    
    
end
legend('model-based','online','data-based')
