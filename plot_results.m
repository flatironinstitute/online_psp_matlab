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


