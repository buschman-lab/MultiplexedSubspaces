function Plot_AutoMotifOnsetDetection(savedir)
%Camden MacDowell - timeless
%see MotifOnset for details. 
%takes ~1 hour to run

if nargin <1; savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\MotifOnset'; end
if~exist(savedir,'dir'); mkdir(savedir); end

fp = fig_params_cortdynamics;

for cur_rec = 1:6
[rec_name,~,~,~,motif_fits] = LoadDataDirectories(cur_rec);

[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%get the motif activity, spatial rho, and H waveform
weight_all = cell(1,numel(motif_fits));
rho_all = cell(1,numel(motif_fits));
H_all = cell(1,numel(motif_fits));
for i = 1:numel(motif_fits)
    data = load(motif_fits{i},'w','H','stats_refit');
    weight = arrayfun(@(n) tensor_convolve(nanmax(data.w(:,n,:),[],1),data.H(n,:)),1:size(data.w,2),'UniformOutput',0);
    weight_all{i} = cat(1,weight{:});
    rho_all{i} = data.stats_refit.rho_frame_per_motif;
    H_all{i} = data.H;
end
weight = [weight_all{:}];
rho = [rho_all{:}];
H = [H_all{:}];

%remove noise motif
weight(fp.noisemotif,:)=[];
rho(fp.noisemotif,:)=[];
H(fp.noisemotif,:)=[];
motif_onset(fp.noisemotif) = [];

%% figures
m = numel(motif_onset);
%H waveform around time of trigger
figure; hold on; 
t=tiledlayout(4,4); t.TileSpacing = 'compact'; t.Padding = 'compact';
for i = 1:m
    nexttile; hold on;
    [~,temp] = ParseByOnset([],{H'},motif_onset,[-5 5],i);
    temp = squeeze(temp{1}(i,:,:)); %just plot this motif
%     temp = temp./max(temp,[],1);
    x = -5:5;    
    shadedErrorBar(x,nanmean(temp'),sem(temp'),'lineprops',{'color',[0.8500 0.3723 0.0078 0.15],'linewidth',2});
%     plot(x,temp,'color',[0.8500 0.3723 0.0078 0.15],'linewidth',0.5);
%     plot(x,nanmean(temp,2),'color','k','linewidth',2);
    yval = get(gca,'ylim');
    plot([0,0],yval,'color','r','linestyle',':','linewidth',2)
    ylim(yval);
    fp.FormatAxes(gca);   
    xvals = get(gca,'xtick');
    set(gca,'XTickLabel',round(xvals*(1000/15)));
    xlabel('time relative to onset (ms)');
    ylabel('H weighting')
    title(sprintf('Motif %d',i),'FontWeight','normal');    
end
title(t,sprintf('%s',rec_name),'FontWeight','normal');
fp.FigureSizing(gcf,[3 2 15 15],[2 2 30 20])

%Spatial correlation around trigger
figure; hold on; 
t=tiledlayout(4,4); t.TileSpacing = 'compact'; t.Padding = 'compact';
for i = 1:m
    nexttile; hold on;
    [~,temp] = ParseByOnset([],{rho'},motif_onset,[-2 20],i);
    x = -2:20; 
    y = squeeze(temp{1}(i,:,:)); %just plot this motif 
    shadedErrorBar(x,nanmean(y'),sem(y'),'lineprops',{'color',[0.8500 0.3723 0.0078 0.55],'linewidth',2});
    y = arrayfun(@(n) squeeze(temp{1}(n,:,:)), find(~ismember(1:m,i)==1),'UniformOutput',0); %plot the others
    y = cat(2,y{:});
    shadedErrorBar(x,nanmean(y'),sem(y'),'lineprops',{'color',[0.75 0.75 0.75 1],'linewidth',0.5});
    yval = get(gca,'ylim');
    plot([0,0],yval,'color','r','linestyle',':','linewidth',2)            
    plot([x(1),x(end)],[0 0],'color','k','linestyle',':','linewidth',2)            
    ylim(yval);
    xlim([x(1),x(end)])
    fp.FormatAxes(gca);   
    xvals = get(gca,'xtick');
    set(gca,'XTickLabel',round(xvals*(1000/15)));
    xlabel('time relative to onset (ms)');
    ylabel('Spatial Correlation')
    title(sprintf('Motif %d',i),'FontWeight','normal');     
end
title(t,sprintf('%s',rec_name),'FontWeight','normal');
fp.FigureSizing(gcf,[3 2 15 15],[2 2 30 20])

%Spatial correlation around trigger
figure; hold on; 
t=tiledlayout(4,4); t.TileSpacing = 'compact'; t.Padding = 'compact';
for i = 1:m
    nexttile; hold on;
    [~,temp] = ParseByOnset([],{rho'},motif_onset,[-2 20],i);
    x = -2:20; 
    y = squeeze(temp{1}(i,:,:)); %just plot this motif 
    shadedErrorBar(x,nanmean(y'),sem(y'),'lineprops',{'color',[0.8500 0.3723 0.0078 0.55],'linewidth',2});
    y = arrayfun(@(n) squeeze(temp{1}(n,:,:)), find(~ismember(1:m,i)==1),'UniformOutput',0); %plot the others
    y = cat(2,y{:});
    shadedErrorBar(x,nanmean(y'),sem(y'),'lineprops',{'color',[0.75 0.75 0.75 1],'linewidth',0.5});
    yval = get(gca,'ylim');
    plot([0,0],yval,'color','r','linestyle',':','linewidth',2)            
    plot([x(1),x(end)],[0 0],'color','k','linestyle',':','linewidth',2)            
    ylim(yval);
    xlim([x(1),x(end)])
    fp.FormatAxes(gca);   
    xvals = get(gca,'xtick');
    set(gca,'XTickLabel',round(xvals*(1000/15)));
    xlabel('time relative to onset (ms)');
    ylabel('Spatial Correlation')
    title(sprintf('Motif %d',i),'FontWeight','normal');     
end
title(t,sprintf('%s',rec_name),'FontWeight','normal');
fp.FigureSizing(gcf,[3 2 15 15],[2 2 30 20])

%Motif activity around trigger
figure; hold on; 
t=tiledlayout(4,4); t.TileSpacing = 'compact'; t.Padding = 'compact';
for i = 1:m
    nexttile; hold on;
    [~,temp] = ParseByOnset([],{weight'},motif_onset,[-2 20],i);
    x = -2:20; 
    y = squeeze(temp{1}(i,:,:)); %just plot this motif 
    shadedErrorBar(x,nanmean(y'),sem(y'),'lineprops',{'color',[0.8500 0.3723 0.0078 0.55],'linewidth',2});
    yval = get(gca,'ylim');
    plot([0,0],yval,'color','r','linestyle',':','linewidth',2)            
    plot([x(1),x(end)],[0 0],'color','k','linestyle',':','linewidth',2)            
    ylim(yval);
    xlim([x(1),x(end)])
    fp.FormatAxes(gca);   
    xvals = get(gca,'xtick');
    set(gca,'XTickLabel',round(xvals*(1000/15)));
    xlabel('time relative to onset (ms)');
    ylabel('Spatial Correlation')
    title(sprintf('Motif %d',i),'FontWeight','normal');     
end
title(t,sprintf('%s',rec_name),'FontWeight','normal');
fp.FigureSizing(gcf,[3 2 15 15],[2 2 30 20])
        
%save off
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('%s_AdaptiveThresholding',rec_name),savedir,0); close all

end %recording loop
end    










