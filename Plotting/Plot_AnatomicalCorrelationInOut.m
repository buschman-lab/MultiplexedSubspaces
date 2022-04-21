function Plot_AnatomicalCorrelationInOut(fast_load,parse_type)
%Camden - timeless
%leave st_mat and st_depth empty to reload from raw
if nargin <1; fast_load = 1; end
if nargin <2; parse_type = 'general'; end

%load ephys data    
for cur_rec = 1:6
    if fast_load==0
        [~,~,~,EphysPath,~] = LoadDataDirectories(cur_rec);
        [st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe'); 
    else %preprocessed
        [rec_name,~,~,EphysPath] = LoadDataDirectories(cur_rec);
        [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDmotif1.mat'],0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA'}); %first motif has all that you need
        temp = cellfun(@(x) load(x,'st_mat','st_depth'),fn);  
        st_depth = temp.st_depth;
        st_mat = temp.st_mat; 
    end
    %split into areas
    neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
    [area_val, area_label] = ParseByArea(cat(2,st_mat{:})',cat(2,neu_area{:}),parse_type);
    [area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 
    %within correlation
    in = cellfun(@(x) nanmean(corr(x'),'all'),area_val);
    out = NaN(numel(area_val),1);
    idx = 1:numel(area_val);
    for i = 1:numel(area_val)            
        out(i) = nanmean(corr(area_val{i}',cat(1,area_val{~ismember(idx,i)})'),'all');
    end
    
    label{cur_rec} = area_label;
    allin{cur_rec} = in;
    allout{cur_rec} = out;
end

x = cat(1,allin{:});
y = cat(1,allout{:});
area_name = cat(1,label{:});
unique_area = unique(area_name);

%parse by area
xx = arrayfun(@(n) nanmean(x(ismember(area_name,unique_area{n}))),1:numel(unique_area));
yy = arrayfun(@(n) nanmean(y(ismember(area_name,unique_area{n}))),1:numel(unique_area));

%plot in vs out, colored by area
fp = fig_params_cortdynamics;
figure; hold on;  
plot([0 max(cat(1,xx(:),yy(:)))],[0 max(cat(1,xx(:),yy(:)))],'linestyle','-','linewidth',1.5,'color','k')
col = getColorPalet(numel(unique_area));
p = cell(1,numel(unique_area));
for i = 1:numel(unique_area)
   p{i} = plot(xx(i),yy(i),'linestyle','none','marker','o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',fp.markersizebig);
end
legend([p{:}],unique_area,'location','northwest');

xlabel('within-area rho')
ylabel('out-of-area rho');
fp.FormatAxes(gca);


end

















