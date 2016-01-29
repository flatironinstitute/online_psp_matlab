%% display results simulation
files=dir('*.mat')
figure
hold all
legends={};
cm=colormap(hot(25))
d_s=[];
q_s=[];
err_s=[];
errhalf_s=[];
for f=1:numel(files)
   disp(files(f).name)
   load(files(f).name,'d','q','errors')
   legends{f}=files(f).name;
   plot(median(errors,2),'d','Linewidth',2,'Color',cm(f,:))
   err_s=[ err_s errors(end,:)'];
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
files=dir('*.mat')
figure
hold all
legends={};
cm=colormap('lines')
% cm=colormap(hot(20))
colq=[];
for f=1:numel(files)
   disp(files(f).name)
   load(files(f).name,'d','q','errors')
   legends{f}=files(f).name;
   legends{f}=['input dim=' num2str(d)];
   colq(f)=errorbar(q+normrnd(0,q/50,size(q)),nanmedian(errors(end,:)),quantile(errors(end,:),.25),quantile(errors(end,:),.75),'ko','MarkerFaceColor',cm(log2(d),:),'MarkerSize',10) ;
   set(gca,'yscale','log')
   set(gca,'xscale','log')
%    legend(legends{:})   
% input('')
end
[vals,idx]=unique(d_s)
legend(colq(idx),legends{idx}, 'Interpreter', 'none')
xlabel(['output dim = [' num2str(unique(q_s)) ']'], 'Interpreter', 'none')
ylabel('Projection Error')
saveas(gcf,'ProjErrors.jpg')
saveas(gcf,'ProjErrors.fig')
%%
files=dir('*.mat')
figure
hold all
legends={};
cm=colormap('lines')
% cm=colormap(hot(20))
colq=[];
for f=1:numel(files)
   disp(files(f).name)
   load(files(f).name,'d','q','times_')
   legends{f}=files(f).name;
   legends{f}=['input dim=' num2str(d)];
   colq(f)=errorbar(q+normrnd(0,q/50,size(q)),nanmedian(times_(:)*1000),quantile(times_(:)*1000,.25),quantile(times_(:)*1000,.75),'ko','MarkerFaceColor',cm(log2(d),:),'MarkerSize',10) ;
   set(gca,'yscale','log')
   set(gca,'xscale','log')
%    legend(legends{:})   
% input('')
end
[vals,idx]=unique(d_s)
legend(colq(idx),legends{idx}, 'Interpreter', 'none')
xlabel(['output dim = [' num2str(unique(q_s)) ']'], 'Interpreter', 'none')
ylabel('Time per iteration (ms)')
saveas(gcf,'TimeIter.jpg')
saveas(gcf,'TimeIter.fig')
%% ********************************** PLOT FOR GAPS *****************************************
files_to_analize = uipickfiles()
%%
figure
hold all
legends={};
cm1=colormap('hot');
cm2=flipud(colormap('gray'))
colq=[];
for ff=1:numel(files_to_analize)
    disp(files_to_analize{ff})
    load(fullfile(files_to_analize{ff},'n4096_d256_q16.mat'),'options_generator','test_method','errors')
    legends{ff}=['rho = ' num2str(options_generator.rho)];
    if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(20+ff,:);
        shift=0;
    else
        symbol='o';
        colr=cm2(ff*2,:);
        shift=0;
    end    
    errline=nanmedian(errors');    
    idx_not_nan=find(~isnan(errline))
    
%     colq(ff)=errorbar(idx_not_nan+shift,nanmedian(errors(idx_not_nan,:)'),mad(errors(idx_not_nan,:)'/1000,1),['-' symbol],'color',colr,'MarkerFaceColor',colr,'MarkerSize',10) ;
    colq(ff)=plot(idx_not_nan+shift,nanmedian(errors(idx_not_nan,:)'),['-' symbol],'color',colr,'MarkerFaceColor',colr,'MarkerSize',10) ;

end
legend(legends)
xlabel('Samples', 'Interpreter', 'none')
ylabel('Projection error')
% set(gca,'yscale','log')
% set(gca,'xscale','log')
%%
files=uipickfiles
%%
figure
hold all
legends={};
cm1=colormap('hot');
cm2=flipud(colormap('gray'))
colq=[];
for f=1:numel(files)
   disp(files{f})
   load(fullfile(files{f},'n4096_d256_q16.mat'),'options_generator','times_','test_method')
   legends{f}=['method=' num2str(test_method)];
   if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(1,:);
        shift=5;
    else
        symbol='o';
        colr=cm1(30,:);
        shift=0;
   end    
   q=options_generator.rho;
   colq(f)=errorbar(q,nanmedian(times_(:)*1000),mad(times_(:)*1000,.25),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
%    legend(legends{:})   
% input('')
end
[vals,idx]=unique(legends)
legend(colq(idx),legends{idx}, 'Interpreter', 'none')
xlabel('rho', 'Interpreter', 'none')
ylabel('Time per iteration (ms)')
saveas(gcf,'TimeIter.jpg')
saveas(gcf,'TimeIter.fig')