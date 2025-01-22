function vp = CompareViolins(data,fp,varargin)
    opts.col = fp.GenDefaultColorArray(size(data,1));
    opts.xpos = (1:1:size(data,1)); 
    opts.label = arrayfun(@(x) num2str(x),1:size(data,1),'UniformOutput',0);
    opts.divfactor = 1;
    opts.connectline = []; %give a color vector if want to plot
    opts.plotspread=0;
    opts.sidebyside = 0; %flag to make it side by side
    opts.plotaverage = 1;
    opts.distWidth = 0.9;
    
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