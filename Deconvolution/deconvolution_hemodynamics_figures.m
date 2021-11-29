function deconvolution_hemodynamics_figures()
%Camden MacDowell - timeless
statsfn = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\hemocorrected\hemocomparestatsstd.mat';
tracefn = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\hemocorrected\\hemocompare_exampletracesstd.mat';

fp = fig_params_deconvolutionpaper;
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\hemocorrection';
if ~exist(savedir,'dir'); mkdir(savedir); end

%plot example traces of the raw signal
traces = load(tracefn);
close all; 
labels = {'none','lr','glm','fNN'};
col = {fp.c_none,fp.c_lr,fp.c_glm,fp.c_ff};
for i = 1:4
    if i ==1; trace_data = traces.nonepred; 
    elseif i ==2; trace_data = traces.lrpred; 
    elseif i ==3; trace_data = traces.glmpred;        
    elseif i ==4; trace_data = cellfun(@(x) cat(1,zeros(floor((numel(traces.glmpred{1})-numel(x))/2),1),x),traces.fNNpred,'UniformOutput',0); %need to add back the offset from fitting the network)
    end
    figure; hold on; 
    for j = 1:2
       x = (1:size(trace_data{j},1))/15;       
       if j==1
           plot(x,trace_data{j},'color','k','linestyle','--','linewidth',1);
       else
           plot(x,trace_data{j},'color',col{i},'linestyle','-','linewidth',1);
       end
    end
    xlim([15,30]); fp.FormatAxes(gca); grid on
    fp.FigureSizing(gcf,[3 2 10 3],[10 10 15 7]); box on        
    saveCurFigs(gcf,{'-dpng','-dsvg'},sprintf('examplepredtrace_%s',labels{i}),savedir,0); close gcf    
    
    %plot the xcorr in the predicted signals    
    [rho,lag] = xcorr(trace_data{1}-nanmean(trace_data{1}),trace_data{2}-nanmean(trace_data{2}),'normalized');    
    figure; hold on; plot(lag,rho,'linewidth',1.5,'color',col{i}); 
    yvals = get(gca,'ylim');
    line([0 0],yvals,'color','k','linewidth',1.5,'linestyle',':'); line([lag(1),lag(end)],[0 0],'color','k','linewidth',1.5,'linestyle',':');
    xlim([-45 45]); set(gca,'ylim',yvals); set(gca,'xtick',[-45:15:45]);
    ylabel('Rho')
    xlabel('Samples');
    fp.FormatAxes(gca); grid on
    fp.SetTitle(gca,labels{i})
    fp.FigureSizing(gcf,[2 2 3 3],[])     
    saveCurFigs(gcf,{'-dpng','-dsvg'},sprintf('predictedxcorr_%szoom',labels{i}),savedir,0); close gcf     
end

%plot the example neural trace
figure; hold on; 
x = (1:size(traces.sttrue,1))/15;       
plot(x,traces.sttrue,'color','k','linestyle','-','linewidth',1);
xlim([15,30]); fp.FormatAxes(gca); grid on
fp.FigureSizing(gcf,[3 2 10 3],[10 10 15 7]); box on        
saveCurFigs(gcf,{'-dpng','-dsvg'},'exampletruetrace',savedir,0); close gcf     
%%
data = load(statsfn);
close all; 
labels = {'none','lr','glm','fNN'};
col = {fp.c_none,fp.c_lr,fp.c_glm,fp.c_ff};
figure; hold on; 
grp = data.grp;
for i = 1:4    
    if i ==1
        temp = [data.stats_none{:}];
        data_rho = cat(1,temp.rho);
    elseif i ==2
        temp = [data.stats_lr{:}];
        data_rho = cat(1,temp.rho);
    elseif i ==3
        temp = [data.stats_glm{:}];
        data_rho = cat(1,temp.rho);       
    elseif i ==4
        temp = [data.stats_fNN{:}];
        data_rho = cat(1,temp.rho);        
    end    
    idx = rand(1,sum(grp==1))/8;
    x = data_rho(grp==1);
    y = data_rho(grp==2);
    b = bar([i,i+0.4],[nanmean(x),nanmean(y)],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor',[1 1 1]);
    b.CData(1,:) = col{i};
    b.CData(2,:) = col{i};
    plot(i+idx-0.06,x,'.','markersize',10,'color','k')
    plot(i+idx+0.4-0.06,y,'.','markersize',10,'color','k')
    %make pairwise
    for j = 1:sum(grp==1)
        plot([i+idx(j)-0.06,i+idx(j)-0.06+0.4],[x(j),y(j)],'color',[0.5 0.5 0.5 0.5],'linewidth',1);
    end
    errorbar(i,nanmean(x),sem(x,1),'LineWidth',1.5,'Color','k');
    errorbar(i+.4,nanmean(y),sem(y,1),'LineWidth',1.5,'Color','k');
    %get the mean and ci of the difference
    d = nanmean(x-y);
    ci = bootci(1000,@nanmean,x-y);
%     p=signrank(fisherZ(x),fisherZ(y));    
    plot([i,i+0.4],[0.6 0.6],'linewidth',1,'color','k');
    text(i+0.2,0.6,sprintf('d=%0.2g \nci=%0.2g-\n%0.2g',d,ci(1),ci(2)),'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',8);       
end
ylabel('Rho');
fp.FormatAxes(gca); grid on
fp.FigureSizing(gcf,[2 2 8 6],[])  
saveCurFigs(gcf,{'-dpng','-dsvg'},'comparisonwhentrue',savedir,0); close gcf   

%% training data version
statsfn = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\hemocorrected\hemocomparestats_trainstd.mat';
tracefn = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\hemocorrected\\hemocompare_exampletraces_trainstd.mat'; 
fp = fig_params_deconvolutionpaper;
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\hemocorrection_train';
if ~exist(savedir,'dir'); mkdir(savedir); end

%plot example traces of the raw signal
traces = load(tracefn);
close all; 
labels = {'none','lr','glm','fNN'};
col = {fp.c_none,fp.c_lr,fp.c_glm,fp.c_ff};
for i = 1:4
    if i ==1; trace_data = traces.nonepred_train; 
    elseif i ==2; trace_data = traces.lrpred_train; 
    elseif i ==3; trace_data = traces.glmpred_train;
    elseif i ==4; trace_data = cellfun(@(x) cat(1,zeros(floor((numel(traces.glmpred_train{1})-numel(x))/2),1),x),traces.fNNpred_train,'UniformOutput',0); %need to add back the offset from fitting the network)
    end
    figure; hold on; 
    for j = 1:2
       x = (1:size(trace_data{j},1))/15;       
       if j==1
           plot(x,trace_data{j},'color','k','linestyle','--','linewidth',1);
       else
           plot(x,trace_data{j},'color',col{i},'linestyle','-','linewidth',1);
       end
    end
    xlim([1070,1085]); fp.FormatAxes(gca); grid on
    fp.FigureSizing(gcf,[3 2 10 3],[10 10 15 7]); box on        
    saveCurFigs(gcf,{'-dpng','-dsvg'},sprintf('examplepredtrace_%s_train',labels{i}),savedir,0); close gcf    
    
    %plot the xcorr in the predicted signals    
    [rho,lag] = xcorr(trace_data{1}-nanmean(trace_data{1}),trace_data{2}-nanmean(trace_data{2}),'normalized');    
    figure; hold on; plot(lag,rho,'linewidth',1.5,'color',col{i}); 
    yvals = get(gca,'ylim');
    line([0 0],yvals,'color','k','linewidth',1.5,'linestyle',':'); line([lag(1),lag(end)],[0 0],'color','k','linewidth',1.5,'linestyle',':');
    xlim([-45 45]); set(gca,'ylim',yvals); set(gca,'xtick',[-45:15:45]);
    ylabel('Rho')
    xlabel('Samples');
    fp.FormatAxes(gca); grid on
    fp.SetTitle(gca,labels{i})
    fp.FigureSizing(gcf,[2 2 3 3],[])     
    saveCurFigs(gcf,{'-dpng','-dsvg'},sprintf('predictedxcorr_%szoom_train',labels{i}),savedir,0); close gcf     
end

%plot the example neural trace
figure; hold on; 
x = (1:size(traces.sttrue_train,1))/15;       
plot(x,traces.sttrue_train,'color','k','linestyle','-','linewidth',1);
xlim([1070,1085]); fp.FormatAxes(gca); grid on
fp.FigureSizing(gcf,[3 2 10 3],[10 10 15 7]); box on  
saveCurFigs(gcf,{'-dpng','-dsvg'},'exampletruetrace',savedir,0); close gcf     

data = load(statsfn);
close all; 
labels = {'none','lr','glm','fNN'};
col = {fp.c_none,fp.c_lr,fp.c_glm,fp.c_ff};
figure; hold on; 
grp = data.grp;
for i = 1:4    
    if i ==1; data_rho = [data.stats_none_train.rho]; 
    elseif i ==2; data_rho = [data.stats_lr_train.rho]; 
    elseif i ==3; data_rho = [data.stats_glm_train.rho];
    elseif i ==4; data_rho = [data.stats_fNN_train.rho];
    end    
    idx = rand(1,sum(grp==1))/8;
    x = data_rho(grp==1);
    y = data_rho(grp==2);
    b = bar([i,i+0.4],[nanmean(x),nanmean(y)],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor',[1 1 1]);
    b.CData(1,:) = col{i};
    b.CData(2,:) = col{i};
    plot(i+idx-0.06,x,'.','markersize',10,'color','k')
    plot(i+idx+0.4-0.06,y,'.','markersize',10,'color','k')
    %make pairwise
    for j = 1:sum(grp==1)
        plot([i+idx(j)-0.06,i+idx(j)-0.06+0.4],[x(j),y(j)],'color',[0.5 0.5 0.5 0.5],'linewidth',1);
    end
    errorbar(i,nanmean(x),sem(x,2),'LineWidth',1.5,'Color','k');
    errorbar(i+.4,nanmean(y),sem(y,2),'LineWidth',1.5,'Color','k');
    %get the mean and ci of the difference
    d = nanmean(x-y);
    ci = bootci(1000,@nanmean,x-y);
%     p=signrank(fisherZ(x),fisherZ(y));    
    plot([i,i+0.4],[0.6 0.6],'linewidth',1,'color','k');
    text(i+0.2,0.6,sprintf('d=%0.2g \nci=%0.2g-\n%0.2g',d,ci(1),ci(2)),'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',8);       
end
ylabel('Rho');
fp.FormatAxes(gca); grid on
fp.FigureSizing(gcf,[2 2 8 6],[])  
saveCurFigs(gcf,{'-dpng','-dsvg'},'comparisonwhentrue_train',savedir,0); close gcf   
end %function end
