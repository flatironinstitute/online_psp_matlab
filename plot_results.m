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
dirFlags = [files.isdir];
files=files(dirFlags);
files = struct2cell(files);
files = files(1,:);

files=uipickfiles;
%%
figure
hold all
legends={};
cm=colormap('lines')
% cm=colormap(hot(20))
colq=[];
for f=1:numel(files)
   disp(files{f})
   load(files{f},'test_method','d','q','errors')
   legends{f}=files{f};
   
   %legends{f}=['input dim=' num2str(d)];
   legends{f}=[test_method ' rho = ' num2str(options_generator.rho) ' d = ' num2str(d) ' n0 = ' num2str(n0)];

   if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(20+ff,:);
        shift=0;
    elseif isequal(test_method,'SGA')
        symbol='o';
        colr=cm2(ff*2,:);
        shift=0;
    else
        symbol='s';
        colr=cm3(ff*2,:);
        shift=0;
    end
   colq(f)=errorbar(q+normrnd(0,q/50,size(q)),nanmedian(errors(end,:)),quantile(errors(end,:),.25),quantile(errors(end,:),.75),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
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
figure
hold all
legends={};
cm=colormap('lines')
% cm=colormap(hot(20))
colq=[];
for f=1:numel(files)
   disp(files{f})
   load(files{f},'test_method','d','q','times_')
   legends{f}=files{f};
   
   %legends{f}=['input dim=' num2str(d)];
   legends{f}=[test_method ' rho = ' num2str(options_generator.rho) ' d = ' num2str(d) ' n0 = ' num2str(n0)];

   if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(20+ff,:);
        shift=0;
    elseif isequal(test_method,'SGA')
        symbol='o';
        colr=cm2(ff*2,:);
        shift=0;
    else
        symbol='s';
        colr=cm3(ff*2,:);
        shift=0;
    end
   colq(f)=errorbar(q+normrnd(0,q/50,size(q)),nanmedian(times_(:)*1000),quantile(times_(:)*1000,.25),quantile(times_(:)*1000,.75),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
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
figure
hold all
legends={};
cm1=colormap('hot');
cm2=flipud(colormap('gray'))
cm3=flipud(colormap('autumn'))

colq=[];
for ff=1:numel(files_to_analize)
    disp(files_to_analize{ff})
    if isdir(files_to_analize{ff})
        load(fullfile(files_to_analize{ff},'n4096_d256_q16.mat'),'options_generator','test_method','errors')
    else
        load(files_to_analize{ff},'options_generator','test_method','errors','d','q','n0')
    end
    legends{ff}=[test_method ' rho = ' num2str(options_generator.rho) ' d = ' num2str(d) ' q = ' num2str(q) ' n0 = ' num2str(n0)];
    if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(20+ff,:);
        shift=0;
    elseif isequal(test_method,'SGA')
        symbol='o';
        colr=cm2(ff*2,:);
        shift=0;
    else
        symbol='s';
        colr=cm3(ff*2,:);
        shift=0;
    end
    errline=nanmedian(errors');    
    idx_not_nan=find(~isnan(errline))
    
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
figure
hold all
legends={};
colq=[];
for f=1:numel(files_to_analize)
   disp(files_to_analize{f})
   if isdir(files_to_analize{ff})
     load(fullfile(files_to_analize{f},'n4096_d256_q16.mat'),'options_generator','times_','test_method','q','d')
   else
     load(files_to_analize{f},'options_generator','times_','test_method','q','d')  
   end
   legends{f}=['method=' num2str(test_method)];
   if isequal(test_method,'IPCA')
        symbol='d';
        colr=cm1(1,:);
        shift=1;
   elseif  isequal(test_method,'SGA')
        symbol='o';
        colr=cm1(30,:);
        shift=1;
   else
        symbol='s';
        colr=cm3(ff*2,:);
        shift=0;
   end    
   
   xvarname='rho';
   xvar=options_generator.rho;
   xvarname='q';
   xvar=q;
   %colq(f)=errorbar(xvar,nanmedian(times_(:)*1000),mad(times_(:)*1000,.25),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;
   colq(f)=errorbar(xvar+ normrnd(0,xvar/20),nanmedian(times_(:)*1000),quantile(times_(:)*1000,.25),quantile(times_(:)*1000,.75),'ko','MarkerFaceColor',colr,'MarkerSize',10) ;

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