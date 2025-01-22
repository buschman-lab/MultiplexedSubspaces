function Figure1j(x,col,area_label)
%x is the local/shared ratio
fp = fig_params_cortdynamics;

%violin plot the ratio
x = reshape(x,size(x,1)*size(x,2),size(x,3));
x = pairedBootstrap(x,@nanmean);

figure; hold on; 
[~,idx] = sort(nanmean(x),'descend');
plot([0,9],[2,2],'color',[0.8 0 0],'LineStyle','--','LineWidth',1)
plot([0,9],[1,1],'color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1)
CompareViolins(x(:,idx)',fp,'col',col(idx),'plotspread',0,'connectline',[],'plotspread',0,'divfactor',0.5,'distWidth',0.85);
set(gca,'Xtick',1:8,'XTickLabel',area_label(idx),'XTickLabelRotation',90);
ylabel('D_{local}/D_{shared}')
ylim([1.4 3.6])  
fp.FormatAxes(gca);
box on; grid off
fp.FigureSizing(gcf,[3 2 7.5 3.75],[10 10 14 10]) 

end