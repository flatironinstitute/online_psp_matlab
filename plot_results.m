%% display results simulation
clear all
files=dir('*.mat')
figure
hold all
legends={};
cm=colormap(hot(numel(files)));
d_s=[];
q_s=[];
err_s=[];
errhalf_s=[];
for f=1:numel(files)
   disp(files(f).name)
   load(files(f).name,'d','q','errors_real')
   legends{f}=files(f).name;
   plot((median(errors_real,2)),'d','Linewidth',2,'Color',cm(f,:))
   err_s=[ err_s errors_real(end,:)'];
   d_s=[d_s d];
   q_s=[q_s q];
   legend(legends{:},'Interpreter', 'none')   
%    input('')
%    cla
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
%%
figure
clear all
files=dir('*.mat')
% dirFlags = [files.isdir];
% files=files(dirFlags);
files = struct2cell(files);
files = files(1,:);
hold all
legends={};
cm1=hot(20+numel(files));
cm2=flipud(gray(numel(files)))
cm3=flipud(autumn(numel(files)))

% cm=colormap(hot(20))
colq=[];
for f=1:numel(files)    
   disp(f)
   load(files{f},'options_algorithm','options_generator','options_simulations','d','q','errors_online')
   errors=errors_online;
   legends{f}=files{f};
   test_method=options_algorithm.pca_algorithm;
   %legends{f}=['input dim=' num2str(d)];
   legends{f}=[test_method ' rho = ' num2str(options_generator.rho) ' d = ' num2str(d) ' n0 = ' num2str(options_simulations.n0)];

   if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(20+f,:);
        colr=cm1(20,:);
        shift=0;
    elseif isequal(test_method,'H_AH_NN_PCA')
        symbol='o';
        colr=cm2(f,:);
        colr=cm2(20,:);
        shift=0;
    else
        symbol='s';
        colr=cm3(f,:);
        shift=0;
   end
   vr= options_generator.rho;
%    colq(f)=errorbar(vr+normrnd(0,vr/50,size(vr)),nanmedian(errors(end,:)),quantile(errors(end,:),.25),quantile(errors(end,:),.75),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
   colq(f)=scatter(vr+normrnd(0,vr/50,size(vr)),nanmedian(errors(end,:)),'ko','MarkerFaceColor',colr) ;
   set(gca,'yscale','log')
   set(gca,'xscale','log')
%    legend(legends{:})   
% input('')
end
[vals,idx]=unique(legends)
legend(colq(idx),legends{idx}, 'Interpreter', 'none')
xlabel(['output dim'], 'Interpreter', 'none')
ylabel('Projection Error')
saveas(gcf,'ProjErrors.jpg')
saveas(gcf,'ProjErrors.fig')
%%
clear all
files=dir('*.mat')
% dirFlags = [files.isdir];
% files=files(dirFlags);
files = struct2cell(files);
files = files(1,:);
cm1=hot(20+numel(files));
cm2=flipud(gray(numel(files)))
cm3=flipud(autumn(numel(files)))
figure
hold all
legends={};
% cm=colormap(hot(20))
colq=[];
for f=1:numel(files)
   disp(f)
   load(files{f},'options_algorithm','options_generator','options_simulations','d','q','times_')
   legends{f}=files{f};
   times_=diff(times_);
   test_method=options_algorithm.pca_algorithm;
   %legends{f}=['input dim=' num2str(d)];
   legends{f}=[test_method ' rho = ' num2str(options_generator.rho) ' d = ' num2str(d) ' n0 = ' num2str(options_simulations.n0)];

   if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(20+f,:);
        shift=0;
    elseif isequal(test_method,'H_AH_NN_PCA')
        symbol='o';
        colr=cm2(f,:);
        shift=0;
    else
        symbol='s';
        colr=cm3(f,:);
        shift=0;
   end
   vr=q;
   colq(f)=errorbar(vr+normrnd(0,vr/50,size(vr)),nanmedian(times_(:)*1000),quantile(times_(:)*1000,.25),quantile(times_(:)*1000,.75),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
   set(gca,'yscale','log')
   set(gca,'xscale','log')
%    legend(legends{:})   
% input('')
end
[vals,idx]=unique(legends);
legend(colq(idx),legends{idx}, 'Interpreter', 'none')
xlabel('output dim' , 'Interpreter', 'none')
ylabel('Time per iteration (ms)')
saveas(gcf,'TimeIter.jpg')
saveas(gcf,'TimeIter.fig')
%% ********************************** PLOT FOR GAPS *****************************************
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
%%
figure
hold all
legends={};
cm1=hot(20+numel(files));
cm2=flipud(gray(numel(files)))
cm3=flipud(autumn(numel(files)))

colq=[];
for ff=1:numel(files_to_analize)
    disp(files_to_analize{ff})   
    load(files{ff},'options_algorithm','options_generator','options_simulations','d','q','errors_online')
    test_method=options_algorithm.pca_algorithm;
    errors=errors_online;
    legends{ff}=[test_method ' rho = ' num2str(options_generator.rho) ' d = ' num2str(d) ' n0 = ' num2str(options_simulations.n0)];
    if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(20+ff,:);
        shift=0;
    elseif isequal(test_method,'H_AH_NN_PCA')
        symbol='o';
        colr=cm2(ff,:);
        shift=0;
    else
        symbol='s';
        colr=cm3(ff,:);
        shift=0;
    end
    errline=nanmedian(errors');    
    idx_not_nan=find(~isnan(errline));
    
%     colq(ff)=errorbar(idx_not_nan+shift,nanmedian(errors(idx_not_nan,:)'),mad(errors(idx_not_nan,:)'/1000,1),['-' symbol],'color',colr,'MarkerFaceColor',colr,'MarkerSize',10) ;
    colq(ff)=plot(idx_not_nan+shift,nanmedian(errors(idx_not_nan,:)'),['-' symbol],'color',colr,'MarkerFaceColor',colr,'MarkerSize',10) ;

end
legend(legends,'Interpreter', 'none')
xlabel('Samples', 'Interpreter', 'none')
ylabel('Projection error')
% set(gca,'yscale','log')
% set(gca,'xscale','log')
saveas(gcf,'Error.jpg')
saveas(gcf,'Error.fig')

%%
clear all
clear all
files=dir('*.mat')
% dirFlags = [files.isdir];
% files=files(dirFlags);
files = struct2cell(files);
files = files(1,:);
cm1=hot(20+numel(files));
cm2=flipud(gray(numel(files)));
cm3=flipud(autumn(numel(files)));
figure
hold all
legends={};
colq=[];
for f=1:numel(files)
   disp(f)  
   load(files{f},'options_algorithm','options_generator','options_simulations','d','q','times_','errors_online')
   legends{f}=files{f};
   times_=diff(times_);
   test_method=options_algorithm.pca_algorithm;
   errors=errors_online;
   if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(100,:);
        shift=1;
   elseif  isequal(test_method,'H_AH_NN_PCA')
        symbol='o';
        colr=cm2(100,:);
        shift=1;
   else
        symbol='s';
        colr=cm3(ff*2,:);
        shift=0;
   end    
   
   xvarname='rho';
   xvar=options_generator.rho;
   xvarname='q';
   xvar=d;
   %colq(f)=errorbar(xvar,nanmedian(times_(:)*1000),mad(times_(:)*1000,.25),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
   times_iter=times_(:)*1000;
%    colq(f)=errorbar(xvar+ normrnd(0,xvar/20),nanmedian(times_iter),quantile(times_iter,.25),quantile(times_iter,.75),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
   colq(f)=errorbar(xvar+ normrnd(0,xvar/20),nanmedian(times_iter),mad(times_iter,1),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;

%    legend(legends{:})   
% input('')
end
[vals,idx]=unique(legends)
legend(colq(idx),legends{idx}, 'Interpreter', 'none')
xlabel(xvarname, 'Interpreter', 'none')
ylabel('Time per iteration (ms)')
set(gca,'yscale','log')
set(gca,'xscale','log')
saveas(gcf,'TimeIter.jpg')
saveas(gcf,'TimeIter.fig')

%%
clear all
clear all
files=dir('*.mat')
files = struct2cell(files);
files = files(1,:);
mu_time=[];
mad_time=[];

mu_err_real=[];
mu_err_batch=[];
mu_err_online=[];
mad_err_real=[];
mad_err_batch=[];
mad_err_online=[];
d_s=[];
q_s=[];
rho_s=[];
methods_={};
for f=1:numel(files)
   disp(f)  
   load(files{f},'options_algorithm','options_generator','options_simulations','d','q','times_','errors_real','errors_batch_pca','errors_online')
   test_method=options_algorithm.pca_algorithm;   
   times_=diff(times_);
   times_iter=times_(:)*1000; 
   
   mu_time(f)=nanmedian(times_iter);
   mad_time(f)=mad(times_iter,1);
   idx_not_nan=find(~isnan(nanmedian(errors_real,2)));
   mu_err_real(f)=nanmedian(errors_real(idx_not_nan(end),:),2);
   mu_err_batch(f)=nanmedian(errors_batch_pca(idx_not_nan(end),:),2);
   mu_err_online(f)=nanmedian(errors_online(idx_not_nan(end),:),2);
   mad_err_real(f)=mad(errors_real(idx_not_nan(end),:),1,2);
   mad_err_batch(f)=mad(errors_batch_pca(idx_not_nan(end),:),1,2);
   mad_err_online(f)=mad(errors_online(idx_not_nan(end),:),1,2);
   d_s(f)=d;
   q_s(f)=q;
   rho_s(f)=options_generator.rho;
   methods_{f}=test_method;      
end
%%
figure
hold all
ax1=[];
ax2=[];
counter=0;
legends={};
cm1=hot(8);
cm2=(gray(10));
xvar=rho_s;
xvar_name='rho';
col_var_name='q';
col_var=q_s;
stats_={'nanmedian',@(X) mad(X,1)};
error=mu_err_real;
for q=unique(col_var)    
    counter=counter+1;    
    legends{counter}=[col_var_name '=' num2str(q)];
    xax=unique(xvar);
    idx=find(col_var==q & strcmp(methods_,'H_AH_NN_PCA'));
    
    [me_h,ma_h]=grpstats(error(idx),xvar(idx),stats_);
    ax1(counter)=errorbar(xax+normrnd(0,.01,size(xax)), me_h,ma_h,'o','MarkerSize',20,'MarkerFaceColor',cm1(counter,:),'color',cm1(counter,:));
    idx=find(col_var==q & strcmp(methods_,'IPCA'));
    [me_i,ma_i]=grpstats(error(idx),xvar(idx),stats_);    
    ax2(counter)=errorbar(xax+normrnd(0,.01,size(xax)), me_i,ma_i,'o','MarkerSize',20,'MarkerFaceColor',cm2(2+counter,:),'color',cm2(2+counter,:));
    
end
legend([ax1 ax2],[legends legends])
xlabel(xvar_name, 'Interpreter', 'none')
ylabel('Projection Error')
%%
figure
q=16;
error=mu_err_real;
stats_={'nanmedian',@(X) mad(X,1)};
idx=find(col_var==q & strcmp(methods_,'H_AH_NN_PCA'));
[me_h,ma_h]=grpstats(error(idx),xvar(idx),stats_);
errorbar(xax+normrnd(0,.01,size(xax)), me_h,ma_h,'o','MarkerSize',20,'MarkerFaceColor',cm1(5,:),'color',cm1(5,:));
hold on
idx=find(col_var==q & strcmp(methods_,'IPCA'));
[me_i,ma_i]=grpstats(error(idx),xvar(idx),stats_);
errorbar(xax+normrnd(0,.01,size(xax)), me_i,ma_i,'o','MarkerSize',20,'MarkerFaceColor',cm2(5,:),'color',cm2(5,:));
legend('H_AH_NN_PCA','IPCA')


