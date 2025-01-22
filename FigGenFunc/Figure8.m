function Figure8(theta,thetacon)

%%
y = {[40, 70], [69,90], [40,90], [40, 70]};
figure;
ViolinAnglePlot(theta{1}(2,:),[0.8 0.1 0.1],1);
ViolinAnglePlot(thetacon{1}(2,:),[0.1 0.1 0.8],2);
set(gca,'xlim',[0.5 2.5])
set(gca,'xtick',[1,2]);
title('WHS - RSP','fontweight','normal')

%%
n = {'WHS-SS','SS-WHS','RSP-WHS'};
for i = 1:3
    figure;
    ViolinAnglePlot(theta{i}(1,:),[0.8 0.1 0.1],1);
    ViolinAnglePlot(thetacon{i}(1,:),[0.1 0.1 0.8],2);
    set(gca,'xlim',[0.5 2.5])
    title(n{i},'fontweight','normal')
    ylim(y{i+1}); 
end



end


function [mu,err] = ViolinAnglePlot(theta,col,xpos)
fp = fig_params_cortdynamics;
CompareViolins(theta,fp,'plotspread',0,'divfactor',1.1,'xpos',xpos,'plotaverage',0,'col',col);
mu = rad2deg(circ_mean(deg2rad(theta)));
line([xpos-0.075 xpos+0.075],[mu mu],'Color',col,'LineWidth',3)
fp.FormatAxes(gca)
xlabel('motif')
fp.FigureSizing(gcf,[3 3 2.5 4],[10 10 14 10])
ylabel('Angle (degrees)');
set(gca,'XTickLabel',{'A','B'},'xlim',[0.5 2.5])
err = rad2deg(circ_std(deg2rad(theta),[],[],2));
hold on; 


end


function vp = CompareViolins(data,fp,varargin)
    opts.col = fp.GenDefaultColorArray(size(data,1));
    opts.xpos = (1:1:size(data,1)); 
    opts.label = arrayfun(@(x) num2str(x),1:size(data,1),'UniformOutput',0);
    opts.divfactor = 1;
    opts.connectline = []; %give a color vector if want to plot
    opts.plotspread=0;
    opts.sidebyside = 0; %flag to make it side by side
    opts.plotaverage = 1;
    opts.distWidth = 0.5;
    
    opts = ParseOptionalInputs(opts,varargin); 
    
    %optionally add pairwise line in background
    if ~isempty(opts.connectline)
        plot(1:size(data,1),nanmedian(data,2),'color',opts.connectline,'linewidth',1.5,'linestyle','-')
    end
    
    %Make plot
    if opts.sidebyside ==1 %only works if two columns        
        vp = distributionPlot(data(1,:)','histOri','left','distWidth',fp.vp_dist_w,'color',opts.col,...
            'histOpt',1.1,'divfactor',opts.divfactor,'addSpread',opts.plotspread,'showMM',0,...
            'xNames',opts.label,'xValues',opts.xpos,'widthDiv',[2 1]);
        hold on
        set(vp{1},'FaceAlpha',fp.vp_alpha);
        vp = distributionPlot(data(2,:)','histOri','right','distWidth',fp.vp_dist_w,'color',opts.col,...
            'histOpt',1.1,'divfactor',opts.divfactor,'addSpread',opts.plotspread,'showMM',0,...
            'xNames',opts.label,'xValues',opts.xpos,'widthDiv',[2 2]);
        hold on
        set(vp{1},'FaceAlpha',fp.vp_alpha);   
        yvals = get(gca,'ylim');
        plot([opts.xpos,opts.xpos],[0 max(yvals)],'linewidth',1,'color',fp.default_color)
   
    elseif opts.sidebyside ==2
        vp = distributionPlot(data','histOri','left','distWidth',fp.vp_dist_w,'color',opts.col,...
        'histOpt',1.1,'divfactor',opts.divfactor,'addSpread',opts.plotspread,'showMM',0,...
        'xNames',opts.label,'xValues',opts.xpos,'distWidth',opts.distWidth);
        hold on
        set(vp{1}(~isnan(vp{1})),'FaceAlpha',fp.vp_alpha);        
        
        if opts.plotaverage
            for i = 1:size(data,1)
                line([opts.xpos(i)-0.075 opts.xpos(i)+0.075],[nanmedian(data(i,:)),nanmedian(data(i,:))],'Color',opts.col{i},'LineWidth',3)
            end
        end
        xlim([0.5 size(data,1)+0.5]);  
    else
        vp = distributionPlot(data','distWidth',fp.vp_dist_w,'color',opts.col,...
            'histOpt',1.1,'divfactor',opts.divfactor,'addSpread',opts.plotspread,'showMM',0,...
            'xNames',opts.label,'xValues',opts.xpos,'distWidth',opts.distWidth);
        hold on
        set(vp{1}(~isnan(vp{1})),'FaceAlpha',fp.vp_alpha);        
        
        if opts.plotaverage
            for i = 1:size(data,1)
                line([opts.xpos(i)-0.075 opts.xpos(i)+0.075],[nanmedian(data(i,:)),nanmedian(data(i,:))],'Color',opts.col{i},'LineWidth',3)
            end
        end
        xlim([0.5 size(data,1)+0.5]);    
    end
    
end