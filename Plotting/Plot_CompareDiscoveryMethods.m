function Plot_CompareDiscoveryMethods(other_methods,cnmf)

%general plot params
xvals = [0,100]; 
yvals = [50,100];
fp = fig_params; 

other_methods = squeeze(struct2cell(other_methods));

%organize data for plotting
num_components = cellfun(@(x) cat(1,x(:).numFactors),other_methods(1,:),'UniformOutput',0);
num_components = cell2mat(num_components{1});

%nmf
spatial_nmf = cellfun(@(x) cat(1,x(:).ExpVar_all),other_methods(1,:),'UniformOutput',0);
spatial_nmf = cellfun(@(x) x*100, spatial_nmf,'UniformOutput',0);

%pca
spatial_pca = cellfun(@(x) x.ExpVar_all,other_methods(2,:),'UniformOutput',0);
spatial_pca = cellfun(@(x) x(:,num_components)*100',spatial_pca,'UniformOutput',0);

%grouping variable
group=[repmat({'nmf'},1,numel(spatial_nmf)),repmat({'pca'},1,numel(spatial_pca))];

%cnmf
cnmf_ev = cat(1,cnmf(:).ExpVar_all)*100; 
cnmf_num = cell2mat(cat(1,cnmf(:).numFactors)); 

%% %plot main 
custom_statfun = @(y)([nanmedian(y);bootci(1000,{@(my)nanmedian(my),y})]); 
clear g
figure('units','centimeters','Position',[5 5 15 15]);
g(1,1)=gramm('x',num_components,'y',[spatial_nmf,spatial_pca],'color',group);
g(1,1).set_layout_options('Position',[0.1 0.1 0.5 0.6],...
    'legend_pos',[0.6 0.64 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);

g(1,1).stat_summary('type',custom_statfun)
g(1,1).axe_property('Xlim',xvals,'Ylim',yvals,'XGrid','on','YGrid','on','YTick',linspace(yvals(1),yvals(2),6));
g(1,1).set_names('x','Dimensions','y','Percent Explained Variance','color','')

%Create cnmf histogram on top
g(2,1)=gramm('x',x,'color',y);
g(2,1).set_layout_options('Position',[0.1 0.7 0.5 0.075],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.1 0.02],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(2,1).stat_density();
g(2,1).axe_property('Xlim',xvals,'XTickLabel','','YTicklabel','','XGrid','on','YGrid','off','TickLength',[0 0]);
g(2,1).set_names('x','','y','')

%Create cnmf histogram on side
g(3,1)=gramm('x',cnmf_ev);
g(3,1).set_layout_options('Position',[0.6 0.1 0.075 0.6],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).stat_density();
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','','YTickLabel','','xlim',yvals,'XGrid','on','YGrid','off','TickLength',[0 0],'XTick',linspace(yvals(1),yvals(2),6));
g(3,1).set_names('x','','y','')

%Set global axe properties
%g.axe_property('color',[0.9 0.9 0.9],'XGrid','on','YGrid','on','GridColor',[1 1 1],'GridAlpha',0.8,'TickLength',[0 0],'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3])
g.axe_property('TickDir','in','GridColor',[0.5 0.5 0.5],'Linewidth',1.5,'FontSize',16,'FontName','Arial','FontWeight','normal');
g.set_color_options('map','brewer_dark');
g.draw();
%%
%reset colors and add line
set(g(3,1).facet_axes_handles.Children,'Color',fp.c_discovery,'LineWidth',2)
set(g(2,1).facet_axes_handles.Children,'Color',fp.c_discovery,'LineWidth',2)
set(g(2,1).facet_axes_handles.YLabel,'string','')
set(g(3,1).facet_axes_handles.YLabel,'string','')
set(g(2,1).facet_axes_handles.Title,'string',{'Motifs Parsimoniously Capture Majority';'of Variance in Neural Activity'},'FontWeight','normal','FontSize',18,'FontName','Arial');
plot(g(1,1).facet_axes_handles,xvals,[nanmedian(cnmf_ev),nanmedian(cnmf_ev)],'Color',fp.c_discovery,'LineWidth',2,'linestyle','--')
plot(g(1,1).facet_axes_handles,[nanmedian(cnmf_num),nanmedian(cnmf_num)],yvals,'Color',fp.c_discovery,'LineWidth',2,'linestyle','--')

