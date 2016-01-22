%% display results simulation
files=dir('*.mat')
close all
hold all
legends={};
cm=colormap(hot(20))
for f=1:numel(files)
   disp(files(f).name)
   load(files(f).name,'d','q','errors')
   legends{f}=files(f).name;
   plot(median(errors,2),'d','Linewidth',2,'Color',cm(f,:))
   legend(legends{:})   
%    input('')
end
legend(legends, 'Interpreter', 'none')
xlabel('Samples')
ylabel('Eigenspace Estimation Error pca vs batch')
%%
files=dir('*.mat')
close all
hold all
legends={};
cm=colormap('lines')
% cm=colormap(hot(20))
colq=[];
for f=1:numel(files)
   disp(files(f).name)
   load(files(f).name,'d','q','times_')
   legends{f}=files(f).name;
   legends{f}=num2str(d);
   colq(f)=errorbar(q+normrnd(0,.02,size(q)),nanmedian(times_(:))*1000,iqr(times_(:)*1000),'ko','MarkerFaceColor',cm(log2(d),:)) ;
   set(gca,'yscale','log')
   set(gca,'xscale','log')
%    legend(legends{:})   
%    input('')
end
legend(colq,legends{:})
xlabel('q')
ylabel('Time per iteration (ms)')

%%
ax=gca;
ax.XTick=[1:numel(legends)];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;


