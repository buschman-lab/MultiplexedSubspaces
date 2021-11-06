function ImpactOfFluorescenceNormalization_GroupwiseStats()
%Camden MacDowell - timeless
%plots the groupwise stats from ImpactOfFluorescenceNormalization
fp = fig_params_deconvolutionpaper;
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Normalization\';
fn = GrabFiles('stats.mat',0,{savedir});
data = cellfun(@(x) load(x),fn,'UniformOutput',0);

%compare the dff vs dff f_stats between hemispheres
dff = cellfun(@(x) x.fstat_between_dff,data,'UniformOutput',1); %dff
dfs = cellfun(@(x) x.fstat_between_dfs,data,'UniformOutput',1); %dfs

%combine for plotting
figure; hold on;
b = bar(1,nanmean(dff),'FaceColor','flat','FaceAlpha',0.7,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.75 0.75 0.75];
b = bar(2,nanmean(dfs),'FaceColor','flat','FaceAlpha',0.7,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.75 0.75 0.75];
%plot the distribution
x = 0.75+rand(numel(dff),1)/2;
for i = 1:numel(x); plot([x(i),x(i)+1],[dff(i),dfs(i)],'linewidth',0.75,'color',[0.1 0.1 0.1 0.5],'marker','.','markersize',9); end
line([1,2],[10 10],'linewidth',2,'color',[0.1 0.1 0.1]);
pval = signrank(dff,dfs);
text(1.5,10,sprintf('%0.3f',pval),'Rotation',0,'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',14);
ylabel('fstat');
set(gca,'units','centimeters','position',[2 2 2.5 4],'xlim',[0.5 2.5]); fp.FormatAxes(gca); 
saveCurFigs(gcf,'-svg','grouped_fstat_methods',savedir,1); close all

%plot the fstat between craniotomy sites with the dff method
in = cellfun(@(x) nanmean(x.fstat_in_dff),data,'UniformOutput',1); %dff
out = cellfun(@(x) nanmean(x.fstat_out_dff),data,'UniformOutput',1); %dff

figure; hold on;
b = bar(1,nanmean(in),'FaceColor','flat','FaceAlpha',0.7,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.75 0.75 0.75];
b = bar(2,nanmean(out),'FaceColor','flat','FaceAlpha',0.7,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.75 0.75 0.75];

%plot individual points
x = 0.75+rand(numel(in),1)/2;
for i = 1:numel(x); plot([x(i),x(i)+1],[in(i),out(i)],'linewidth',0.75,'color',[0.1 0.1 0.1 0.5],'marker','.','markersize',9); end
line([1,2],[1.5 1.5],'linewidth',2,'color',[0.1 0.1 0.1]);
pval = signrank(in,out);
text(1.5,1.5,sprintf('%0.3f',pval),'Rotation',0,'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',14)
ylabel('fstat');
set(gca,'units','centimeters','position',[2 2 3 4],'xlim',[0.5 2.5]); fp.FormatAxes(gca); 
saveCurFigs(gcf,'-svg','grouped_fstat_betweensites',savedir,1); close all






%Tonight/Tomorrow morning: 
%retrace probe locations (do them perfectly)
%make the figures for each animal
%then write the function for all of the ephys as the final plot
%write function to plot on other peoples available data 


% %Compare Ephys 
% load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\restingstate_processed_fn.mat','spike_opts_list') 
% params.mua = 0; %1= use both 'good' and 'mua' units. 0 = just 'good'
% params.depth = [0 600]; %depth from surface of probe
% params.radius = 2; %pixel radius around probe tip    
% params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
% params.bindata = 0; %1=temporally bin data (i.e. to frame rate 15 from 30)
% [~,st,n_neurons] = CompileData_deconvolution([],spike_opts_list(3),params); 
% 
% %convert to firing rate (hZ) and normalize for the number of neurons and get the std
% fr_std = arrayfun(@(n) nanstd(30*(st{n}./n_neurons(n,:)),[],1),1:numel(st),'UniformOutput',0);
% fr_std = cat(1,fr_std{:});
% 
% grp = [1:4].*ones(numel(st),1);
% anovan(fr_std(:),grp(:));
% 




