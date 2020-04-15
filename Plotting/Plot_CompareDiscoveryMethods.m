function [stats] = Plot_CompareDiscoveryMethods(other_methods,cnmf,stnmf)

%general plot params
xvals = [0,125]; 
yvals = [50,100];
fp = fig_params; 

other_methods = squeeze(struct2cell(other_methods));
%%
%organize data for plotting
num_components = 1:numel(other_methods{1,1});

%stnmf
st_nmf_data = {stnmf(:).pev};
temp_dim = {stnmf(:).dim};

%linearly interpolate to match sample rate of sNMF. 
st_nmf_data = arrayfun(@(n) interp1(temp_dim{n},st_nmf_data{n},1:max(temp_dim{n})),1:numel(st_nmf_data),'UniformOutput',0);
st_nmf_data = cellfun(@(x) x'*100,st_nmf_data,'UniformOutput',0);

%nmf
spatial_nmf = cellfun(@(x) cat(1,x(:).ExpVar_all), other_methods(1,:),'UniformOutput',0);
spatial_nmf = cellfun(@(x) x*100, spatial_nmf,'UniformOutput',0);

%pca
spatial_pca = cellfun(@(x) x.ExpVar_all, other_methods(2,:),'UniformOutput',0);
spatial_pca = cellfun(@(x) (x(:,num_components)*100)',spatial_pca,'UniformOutput',0);

%grouping variable
group = [repmat({'stNMF'},1,numel(st_nmf_data)),repmat({'sPCA'},1,numel(spatial_pca)),repmat({'sNMF'},1,numel(spatial_nmf))];

%cnmf
cnmf_ev = cat(1,cnmf(:).ExpVar_all)*100; 
cnmf_dim = cell2mat(cat(1,cnmf(:).numFactors)); 

% group_sub = [repmat({'sPCA'},1,numel(spatial_pca)),repmat({'sNMF'},1,numel(spatial_nmf)),repmat({'stNMF'},1,numel(st_nmf_pev)),repmat({'Motifs'},1,numel(cnmf_dim))];
group_sub = [repmat({'stNMF'},1,numel(st_nmf_data)),repmat({'sPCA'},1,numel(spatial_pca)),repmat({'sNMF'},1,numel(spatial_nmf)),repmat({'Motifs'},1,numel(cnmf_dim))];

%paired distribution of dimensions at exp var of motifs and 
spatial_pca_dim = NaN(1,numel(cnmf));
spatial_nmf_dim = NaN(1,numel(cnmf));
spatial_pca_ev = NaN(1,numel(cnmf));
spatial_nmf_ev = NaN(1,numel(cnmf));
stnmf_dim=NaN(1,numel(cnmf));
stnmf_ev=NaN(1,numel(cnmf));

for i = 1:numel(cnmf)
    spatial_pca_dim(i) = find(other_methods{2,i}.ExpVar_all>=cnmf(i).ExpVar_all,1,'first');
    spatial_nmf_dim(i) = find([other_methods{1,i}.ExpVar_all]>=cnmf(i).ExpVar_all,1,'first');
    if ~isempty(find(st_nmf_data{i}>=(cnmf(i).ExpVar_all*100),1,'first')) %for a few fits, even allowing 250 fits is not enough. 
        stnmf_ev(i) = st_nmf_data{i}(cell2mat(cnmf(i).numFactors));
        stnmf_dim(i) = find(st_nmf_data{i}>=(cnmf(i).ExpVar_all*100),1,'first');
    else
        stnmf_ev(i) = NaN;
        stnmf_dim(i)=NaN;        
    end
    temp = [other_methods{2,i}.ExpVar_all];
    spatial_pca_ev(i) = temp(cell2mat(cnmf(i).numFactors));
    temp = [other_methods{1,i}.ExpVar_all];
    spatial_nmf_ev(i) = temp(cell2mat(cnmf(i).numFactors));    
end
    
%distribution of exp var at dimensions
% spatial_nmf_ev = cellfun(@(x) x(median(cnmf_dim)), spatial_nmf,'UniformOutput',0);
% spatial_pca_ev = cellfun(@(x) x(median(cnmf_dim)), spatial_pca,'UniformOutput',0);


%% %plot main 
custom_statfun = @(y)([nanmedian(y);bootci(1000,{@(my)nanmedian(my),y})]); 
clear g
figure('units','centimeters','Position',[5 5 15 15]);
g(1,1)=gramm('x',num_components,'y',cat(2,st_nmf_data,spatial_pca,spatial_nmf,{zeros(1,numel(num_components))}),'color',[group,{'Motifs'}]);
g(1,1).set_layout_options('Position',[0.1 0.1 0.5 0.6],...
    'legend_pos',[0.6 0.64 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);

g(1,1).stat_summary('type',custom_statfun)
g(1,1).axe_property('Xlim',xvals,'Ylim',yvals,'XGrid','on','YGrid','on','YTick',linspace(yvals(1),yvals(2),6));
g(1,1).set_names('x','Dimensions','y','Percent Explained Variance','color','')


%Create cnmf histogram on top
g(2,1)=gramm('x',cat(1,stnmf_dim',spatial_pca_dim',spatial_nmf_dim',cnmf_dim),'color',group_sub);
g(2,1).set_layout_options('Position',[0.1 0.7 0.5 0.075],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.1 0.02],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(2,1).stat_density('npoints',500);
g(2,1).axe_property('Xlim',xvals,'XTickLabel','','YTicklabel','','XGrid','on','YGrid','off','TickLength',[0 0]);
g(2,1).set_names('x','','y','')

%Create cnmf histogram on side
g(3,1)=gramm('x',cat(1,stnmf_ev',spatial_pca_ev'*100,spatial_nmf_ev'*100,cnmf_ev),'color',[group_sub]);
g(3,1).set_layout_options('Position',[0.6 0.1 0.075 0.6],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).stat_density('npoints',500);
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','','YTickLabel','','xlim',yvals,'XGrid','on','YGrid','off','TickLength',[0 0],'XTick',linspace(yvals(1),yvals(2),6));
g(3,1).set_names('x','','y','')

%Set global axe properties
%g.axe_property('color',[0.9 0.9 0.9],'XGrid','on','YGrid','on','GridColor',[1 1 1],'GridAlpha',0.8,'TickLength',[0 0],'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3])
g.axe_property('TickDir','in','GridColor',[0.3 0.3 0.3],'Linewidth',1.5,'FontSize',16,'FontName','Arial','FontWeight','normal');
g(1,1).set_color_options('map',[fp.c_discovery;fp.c_nmf;fp.c_pca;fp.c_spacetime]);
g(2,1).set_color_options('map',[fp.c_discovery;fp.c_nmf;fp.c_pca;fp.c_spacetime]);
g(3,1).set_color_options('map',[fp.c_discovery;fp.c_nmf;fp.c_pca;fp.c_spacetime]);
g.draw();
%%
%reset colors and add line
set(g(2,1).facet_axes_handles.YLabel,'string','')
set(g(3,1).facet_axes_handles.YLabel,'string','')
set(g(2,1).facet_axes_handles.Title,'string',{'Motifs Parsimoniously Capture Majority';'of Variance in Neural Activity'},'FontWeight','normal','FontSize',18,'FontName','Arial');
plot(g(1,1).facet_axes_handles,xvals,[nanmedian(cnmf_ev),nanmedian(cnmf_ev)],'Color',fp.c_discovery,'LineWidth',2,'linestyle','--')
plot(g(1,1).facet_axes_handles,[nanmedian(cnmf_dim),nanmedian(cnmf_dim)],yvals,'Color',fp.c_discovery,'LineWidth',2,'linestyle','--')

%% Add significance testing on
stats.median_dim = nanmedian(cnmf_dim);
stats.ci_median_dim = bootci(1000,@nanmedian,cnmf_dim);
stats.median_ev = nanmedian(cnmf_ev);
stats.ci_median_ev = bootci(1000,@nanmedian,cnmf_ev);
[stats.pval_spca_ev,~] = signrank(cnmf_ev,spatial_pca_ev);
[stats.pval_snmf_ev,~] = signrank(cnmf_ev,spatial_nmf_ev);
[stats.pval_stnmf_ev,~] = signrank(cnmf_ev,stnmf_ev);
[stats.pval_spca_dim,~] = signrank(cnmf_dim,spatial_pca_dim);
[stats.pval_snmf_dim,~] = signrank(cnmf_dim,spatial_nmf_dim);
[stats.pval_stnmf_dim,~] = signrank(cnmf_dim,stnmf_dim);
stats.median_spca_dim = nanmedian(spatial_pca_dim);
stats.median_nmf_dim = nanmedian(spatial_nmf_dim);
stats.median_stnmf_dim = nanmedian(stnmf_dim);
stats.median_spca_ev = nanmedian(spatial_pca_ev);
stats.ci_median_spca_ev = bootci(1000,@nanmedian,spatial_pca_ev);
stats.median_snmf_ev = nanmedian(spatial_nmf_ev);
stats.ci_median_snmf_ev = bootci(1000,@nanmedian,spatial_nmf_ev);
stats.median_stnmf_ev = nanmedian(stnmf_ev);
stats.ci_median_stnmf_ev = bootci(1000,@nanmedian,stnmf_ev);








































