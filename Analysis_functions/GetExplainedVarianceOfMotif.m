function GetExplainedVarianceOfMotif()

[fn,~] = GrabFiles('\w*train\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});

data = cellfun(@(x) load(x,'stats_refit'),fn,'UniformOutput',0);
pev = cellfun(@(x) x.stats_refit.pev,data,'UniformOutput',1);

figure; hold on; 
histogram(pev,'BinWidth',.01,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'facealpha',0.8)

nanmean(pev(:))
[~,stat2s] = pairedBootstrap(pev(:),@nanmean)


end
