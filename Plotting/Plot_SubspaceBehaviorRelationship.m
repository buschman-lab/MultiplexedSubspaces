function Plot_SubspaceBehaviorRelationship()
fp = fig_params_cortdynamics;
data_all = LoadBehavioralNetworks(0);

%Loop through recs and grab all motifs and areas
area = cell(1,6);
x_all = cell(1,6);
y_all = cell(1,6);
for cur_rec = 1:6
    data = data_all{cur_rec};

    x = arrayfun(@(n) abs(data(n).rho_all),1:size(data,2),'UniformOutput',0); %whisk, nose, shoulder
    y = arrayfun(@(n) abs(data(n).rho_all),1:size(data,2),'UniformOutput',0); %whisk, nose, shoulder
    xsig = arrayfun(@(n) data(n).sig_thresh,1:size(data,2),'UniformOutput',0); %whisk, nose, shoulder
    %threshold
    for j= 1:numel(y)
       y{j}(y{j}<xsig{j})=NaN;
    end
    area{cur_rec} = arrayfun(@(n) abs(data(n).cur_a),1:size(data,2),'UniformOutput',1);
    x_all{cur_rec} = cat(3,x{:});
    y_all{cur_rec} = cat(3,y{:});
end

label = {'whisker pad','nose','shoulder'};

%plot average weightings across areas
a = cat(2,area{:});
y = cat(3,x_all{:});
for j = 1:3
    figure; hold on; 
    for i = 1:8
        col = fp.c_area;
        temp = squeeze(y(j,:,a==i))';
        shadedErrorBar(1:10,nanmean(temp),sem(temp,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});
    end
    title(sprintf('%s',label{j}),'fontweight','normal');
    legend(data_all{1}(1).area_label,'location','bestoutside')    
    fp.FormatAxes(gca);
    fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])    
end

%plot significant weightings across areas
y = cat(3,y_all{:});
for j = 1:3
    figure; hold on; 
    for i = 1:8
        col = fp.c_area;
        temp = squeeze(y(j,:,a==i))';
        shadedErrorBar(1:10,nanmean(temp),sem(temp,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});
    end
    title(sprintf('%s',label{j}),'fontweight','normal');
    legend(data_all{1}(1).area_label,'location','bestoutside')    
    fp.FormatAxes(gca);
    fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])
end

%% box plot version
%plot average weightings across areas
a = cat(2,area{:});
y = cat(3,y_all{:});
for j = 1:3
    x = squeeze(y(j,:,:))';
    figure; hold on;
    boxplot(x,'Notch','on','color',[0.5 0.5 0.5]);
    b = get(get(gca,'children'),'children');
    t = get(b,'tag');
    idx = find(strcmp(t,'Box')==1);
    for i = 1:numel(idx)
    set(b(idx(i)), 'Color', [0.5 0.5 0.5],'LineWidth',1.5);
    end
    xlabel('Dimension');
    ylabel('Rho');
    for i = 1:8
        col = fp.c_area;
        temp = squeeze(y(j,:,a==i))';
        plot(1:10,nanmean(temp),'color',[col(i,:),0.50],'linewidth',2);
    end
    legend(data_all{1}(1).area_label,'location','bestoutside') 
    fp.FormatAxes(gca);  box on; grid on
    fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 14])
    yyaxis right
    plot(100*sum(isnan(x)/size(x,1)),'linewidth',2,'color','k')
    set(gca,'YColor','k');
    ylabel({'% of subspaces with','significant behavior relationship'})
    [p,~,stats] = kruskalwallis(x,[],'off');
    title(sprintf('Relationship to %s motion energy \np=%0.3d',label{j},p),'fontweight','normal')
    figure; 
    multcompare(stats);
end

end