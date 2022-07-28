function stats = PlotLOOsubspaceFit()

%load the data
[fn,~] = GrabFiles('\w*LOO\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace'});

%remove files corresponding to motif 2 (noise)
fn(cellfun(@(x) ~isempty(regexp(x,'motif2','match')),fn,'UniformOutput',1))=[];

%load all  differences in rsq
d = cellfun(@(x) load(x,'delta_rsq_rel'),fn,'UniformOutput',0); %absolute rsq
d = cellfun(@(x) x.delta_rsq_rel*100,d,'UniformOutput',0);

%fyi - not breaking it out by area yet, because the value doesn't tell us a
%lot (dependent on number of neurons, etc). Provided they all are less than
%full model

%loop through dimensions
data = cell(1,30);
for i = 1:30
   a = cellfun(@(x) x(:,:,i),d,'UniformOutput',0);
   a = cellfun(@(x) x(:), a,'UniformOutput',0);
   %combine across recordings and motifs
   a = cat(1,a{:});
   data{i} = a;
end
data = cat(2,data{:});

fp = fig_params_cortdynamics;
figure; hold on; 
shadedErrorBar(1:30,nanmean(data),sem(data),'lineprops',{'color',[0.25 0.25 0.25 0.75],'linewidth',2},'patchSaturation',0.075);
ylabel({'Percentage of','Explainable Variance'});
xlabel('Subspace Dimensions');
fp.FormatAxes(gca); box on; grid on
fp.FigureSizing(gcf,[3 2 5 5],[2 10 10 10])
title('PEV includes non-contributory (i.e. <=0) models');        

%for each dimensions test the fraction that did not contribute
no_contrib = 100*(sum(data<=0)./sum(~isnan(data)));
figure; hold on; 
posdata = data;
posdata(data<=0)=NaN;
shadedErrorBar(1:30,nanmean(posdata),sem(posdata),'lineprops',{'color',[0.25 0.25 0.25 0.75],'linewidth',2},'patchSaturation',0.075);
ylabel({'Percentage of','Explainable Variance'});
xlabel('Subspace Dimensions');
ylim([0 10]);
yyaxis right
plot(1:30,no_contrib,'k','linewidth',2);
ylabel({'Percent non-','contributory models'});
set(gca,'YColor','k');
xlim([0 30])
ylim([0 20])
fp.FormatAxes(gca); box on; grid on
fp.FigureSizing(gcf,[3 2 5 5],[2 10 10 10])
title({'PEV excludes','non-contributory models'},'fontweight','normal');        

% %get the average contribution +/- SEM at the last dimension
% stats.mean = nanmean(posdata(:,end));
% stats.sem = sem(posdata(:,end));
% stats.ci = bootci(1000,@nanmean,posdata(:,end));
% 
% %Get the full model
% d = cellfun(@(x) load(x,'delta_full'),fn,'UniformOutput',0); %absolute rsq
% d = cellfun(@(x) x.delta_full*100,d,'UniformOutput',0);
% data = cellfun(@(x) x(:),d,'UniformOutput',0);
% data = cat(1,data{:});
% 
% %get the average contribution for those that participated
% posdata = data(data>0);
% stats.meanfull = nanmean(posdata);
% stats.cifull = bootci(1000,@nanmean,posdata);
% stats.nocontribfull = numel(posdata)/sum(~isnan(data));





end