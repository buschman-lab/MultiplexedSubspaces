function AddAreaPSTH(area_val,area_label,offset,offsettxt)
fp = fig_params_cortdynamics;
for i = 1:numel(area_label)    
    if i>1
        temp = cat(1,area_val{1:i-1});
        idx = [size(temp,1)+1,size(temp,1)+size(area_val{i},1)];
    else
        idx = [1,size(area_val{i},1)];
    end
    plot([offset offset],idx,'color',fp.c_area(i,:),'linewidth',3); 
    text(offsettxt, idx(1)+(idx(2)-idx(1))/2, area_label{i},'FontWeight','bold','HorizontalAlignment','right','Color',fp.c_area(i,:),'fontsize',fp.font_size,'FontName',fp.font_name)
    
end

end