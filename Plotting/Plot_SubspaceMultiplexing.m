function Plot_SubspaceMultiplexing(data)
%Camden MacDowell - timeless
%looks trial-to-trial on the relationship between subspace dimensions
fp = fig_params_cortdynamics;

%get example data from recording 1
[rho,~,~,~,rsq] = SubspaceDim_Trial(data{6},0);
[~,area_label] = LoadVariable(data,'rel_performance',[]);

%% plot the average activity 
close all;
for a = 1:8
    x = rsq(a,6); %one area across motifs
    x = cat(1,x{:});
    x = sort(x,2,'descend');
    x =  cumsum(x,2)./sum(x,2);
    col = fp.c_area; col = col(strcmp(area_label,area_label{a}),:);
    figure; hold on;
    arrayfun(@(n) plot(1:size(x,2), x(n,:),'linestyle','-','marker','none','color',[0.85 0.85 0.85],'linewidth',0.5,'markersize',4), 1:size(x,1));
    plot(1:size(x,2), nanmean(x),'linestyle','-','marker','o','color',col,'linewidth',1.5,'markersize',5,'MarkerFaceColor','none');  
    % shadedErrorBar(1:size(x,2),nanmean(x),std(x,1),'lineprops',{'color',[col,0.25],'linewidth',2});
    xlim([1 size(x,2)]); ylim([0 max(x(:))])
    plot([0 size(x,2)],[0.8 0.8] ,'linestyle','--','color','k','linewidth',1.5);
    xlabel('# of dimensions (sorted)');
    ylabel({'Contribution to trial (r^2)'});
    title(sprintf('rec 1\n %s | motif %d',area_label{a},5),'fontweight','normal')
    fp.FormatAxes(gca);   box on; grid on
    fp.FigureSizing(gcf,[3 2 3.5 3.5],[10 10 10 10])
    set(gca,'xtick',[1,5,10,15]);
end

% After the example, plot for each areas
figure; hold on;
for a = 1:8
    x = rsq(a,:); %one area across motifs
    x = cat(1,x{:});
    x = sort(x,2,'descend');
    x =  cumsum(x,2)./sum(x,2);
    col = fp.c_area; col = col(strcmp(area_label,area_label{a}),:);
    shadedErrorBar(1:size(x,2),nanmean(x,1),sem(x,1),'lineprops',{'color',[col 0.75],'linewidth',2});
end
set(gca,'xtick',[1,5,10,15]);
xlim([1 size(x,2)]); ylim([0 max(x(:))])
plot([0 size(x,2)],[0.8 0.8] ,'linestyle','--','color','k','linewidth',1.5);
xlabel('# of dimensions (sorted)');
ylabel({'Contribution to trial (r^2)'});
title(sprintf('rec 1\n %s | motif %d',area_label{a},5),'fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 3.5 3.5],[10 10 10 10])


%% plot stability across recordings and areas
% rho = SubspaceDim_Trial(data{1},0);
% rho = fisherZ(rho);
% x = pairedBootstrap(rho',@nanmean);
% [~,ind] = sort(nanmean(x),'descend');
% figure; hold on;       
% col = fp.c_area;
% col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
% %sort
% col = col(ind);
% CompareViolins(x(:,ind)',fp,'plotspread',0,'divfactor',0.5,'plotaverage',1,'col',col,'distWidth',0.75);
% set(gca,'XTickLabel',area_label(ind),'XTickLabelRotation',45)
% yval = get(gca,'ylim');
% ylabel('% Generalization');
% fp.FormatAxes(gca);  box on; grid on
% fp.FigureSizing(gcf,[3 3 6 3],[10 10 10 10])
%% get the bootstrapped stats for consistency of this pattern across recordings/motifs





