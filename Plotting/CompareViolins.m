function vp = CompareViolins(data,fp,varargin)
    opts.col = fp.GenDefaultColorArray(size(data,1));
    opts.xpos = (1:1:size(data,1)); 
    opts.label = arrayfun(@(x) num2str(x),1:size(data,1),'UniformOutput',0);
    opts.divfactor = 1;
    
    opts = ParseOptionalInputs(opts,varargin); 
    
    %Make plot
    vp = distributionPlot(data','distWidth',fp.vp_dist_w,'color',opts.col,...
        'histOpt',1.1,'divfactor',opts.divfactor,'addSpread',0,'showMM',0,...
        'xNames',opts.label,'xValues',opts.xpos);
    hold on
    set(vp{1},'FaceAlpha',fp.vp_alpha);
    
    for i = 1:size(data,1)
        line([opts.xpos(i)-0.075 opts.xpos(i)+0.075],[nanmedian(data(i,:)),nanmedian(data(i,:))],'Color',opts.col{i},'LineWidth',3)
    end
    xlim([0.5 size(data,1)+0.5]);
end