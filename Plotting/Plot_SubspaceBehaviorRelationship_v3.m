function Plot_SubspaceBehaviorRelationship_v3()
fp = fig_params_cortdynamics;
data_all = LoadBehavioralNetworks();

%Loop through recs and grab all motifs and areas
data = cat(2,data_all{:});
rho = abs(cat(1,data.rho_all));
rho_perm = fisherZ(cat(1,data.rho_perm_all));
lag = cat(1,data.lag_all);

% label = {'whisker pad','nose','shoulder'};

%% Full version violin
temp = fisherZ(rho); 

%get the FDR corrected significance    
p = prctile(rho_perm(:,1),95);
prctile(rho_perm(:),90);

figure; hold on; 
[a,stats] = pairedBootstrap(temp,@nanmean);
CompareViolins(a',fp,'plotspread',0,'divfactor',1,'plotaverage',1,'col',repmat({[0.15 0.15 0.15]},1,10),'distWidth',0.75,'connectline',[0.75 0.75 0.75]);
plot([0 10.5],[p p],'linestyle',':','color',[0.8 0.1 0.1],'linewidth',1.5)
xlim([0.5 10.5])
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 5 4],[2 10 10 10])
xlabel('Subspace dimension');
ylabel('rho_z')
ylim([0.05 0.35])

saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'Correlation_Full','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\BehavioralRelationship',0); close all




%% same for the size
maps = LoadCorticalNetworks(0);
map = cat(2,maps{:});
%go through each one
x = NaN(size(map,2),10);
for j = 1:size(map,2)
    y = map(j).rho_all;
    p = map(j).sig_thresh;    
    p = map(j).sig_thresh;    
    for i = 1:10
        temp = abs(y(:,:,i));
        x(j,i) = 100*(sum(temp>p(i),'all')/sum(~isnan(temp(:))));        
    end
end
% show that the first dimension spans the cortex
figure; hold on; 
[a,stats] = pairedBootstrap(x,@nanmean);
CompareViolins(a',fp,'plotspread',0,'divfactor',1,'plotaverage',1,'col',repmat({[0.15 0.15 0.15]},1,10),'distWidth',0.75,'connectline',[0.75 0.75 0.75]);
xlim([0.5 10.5])
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 5 4],[2 10 10 10])
xlabel('Subspace dimension');
ylabel('Cortical map size (%)')
ylim([15 80])

saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'MapSize','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\BehavioralRelationship',0); close all





%%
% 
% 
% 
% %% Full version bar
% temp = fisherZ(rho); 
% 
% %get the FDR corrected significance    
% p = prctile(rho_perm(:,1),95);
% prctile(rho_perm(:),90);
% 
% figure; hold on; 
% for i = 1:10
%     a = temp(:,i);
%     ci = bootci(1000,@nanmean,a(:)); a = nanmean(a);
%     errorbar(i,a,a-ci(1),a-ci(2),'linestyle','none','marker','o',...
%         'MarkerEdgeColor','none','MarkerFaceColor',[0.5 0.5 0.5],'markersize',5,'linewidth',1.5,'color',[0.5 0.5 0.5],'CapSize',10)
% end
% plot(nanmean(temp),'color',[0.5 0.5 0.5],'linewidth',1);
% plot([0 10.5],[p p],'linestyle',':','color',[0.8 0.1 0.1],'linewidth',1.5)
% xlim([0.5 10.5])
% fp.FormatAxes(gca); box on; grid on; 
% fp.FigureSizing(gcf,[3 2 5 5],[2 10 10 10])
% xlabel('dimensions');
% ylabel('rho_z')
% ylim([0.05 0.35])
% 
% 
% 
% %% Plot non-bootstrapped version
%       
% temp = fisherZ(rho);
% a = temp(:,1);
% b = temp(:,2:3); 
% c = temp(:,4:end);  
% 
% %get the FDR corrected significance    
% prctile(rho_perm(:),95)
% prctile(rho_perm(:),90)
% 
% figure; hold on; 
% ci = bootci(1000,@nanmean,a(:)); a = nanmean(a);
% errorbar(1,a,a-ci(1),a-ci(2),'linestyle','none','marker','o',...
%     'MarkerEdgeColor','none','MarkerFaceColor',[0.5 0.5 0.5],'markersize',6,'linewidth',1.5,'color',[0.5 0.5 0.5],'CapSize',10)
% ci = bootci(1000,@nanmean,b(:)); b = nanmean(b(:));
% errorbar(2,b,b-ci(1),b-ci(2),'linestyle','none','marker','o',...
%     'MarkerEdgeColor','none','MarkerFaceColor',[0.5 0.5 0.5],'markersize',6,'linewidth',1.5,'color',[0.5 0.5 0.5],'CapSize',10)
% ci = bootci(1000,@nanmean,c(:)); c = nanmean(c(:));
% errorbar(3,c,c-ci(1),c-ci(2),'linestyle','none','marker','o',...
%     'MarkerEdgeColor','none','MarkerFaceColor',[0.5 0.5 0.5],'markersize',6,'linewidth',1.5,'color',[0.5 0.5 0.5],'CapSize',10)
% 
% xlim([0.5 3.5])
% fp.FormatAxes(gca); box on; grid on; 
% fp.FigureSizing(gcf,[3 2 3 5],[2 10 10 10])
% set(gca,'xtick',[1,2,3]);
% xlabel('dimensions');
% ylabel('rho_z')
% ylim([0 1])
% % title('Motor Activity','fontweight','normal');
% %     [~,stats] = pairedBootstrap(temp,@nanmean);
% 
% 
% %% bootstrapped version    
% temp = fisherZ(rho);
% a = temp(:,1);
% a = pairedBootstrap(a,@nanmean);
% b = temp(:,2:3);
% b = pairedBootstrap(b(:),@nanmean);  
% c = temp(:,4:end);
% c = pairedBootstrap(c(:),@nanmean);    
% 
% %get the FDR corrected significance    
% prctile(rho_perm(:),95)
% prctile(rho_perm(:),90)
% 
% figure; hold on; 
% CompareViolins(a',fp,'plotspread',0,'divfactor',1,'plotaverage',1,'col',{fp.c_ff},'distWidth',0.75,'xpos',1);
% CompareViolins(b(:)',fp,'plotspread',0,'divfactor',1,'plotaverage',1,'distWidth',0.75,'xpos',2);    
% CompareViolins(c(:)',fp,'plotspread',0,'divfactor',1,'plotaverage',1,'distWidth',0.75,'xpos',3); 
% 
% xlim([0.5 3.5])
% fp.FormatAxes(gca); box on; grid on; 
% fp.FigureSizing(gcf,[3 2 3 5],[2 10 10 10])
% set(gca,'xtick',[1,2,3]);
% xlabel('dimensions');
% ylabel('rho_z')
% ylim([0 1])
% title(label{j},'fontweight','normal');
% %     [~,stats] = pairedBootstrap(temp,@nanmean);
% 
% 
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'Trial_Correlation','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\BehavioralRelationship',0); close all

%% Go ahead and redo the size map as well in the same format
% 
% maps = LoadCorticalNetworks(0);
% map = cat(2,maps{:});
% %go through each one
% x = NaN(size(map,2),10);
% for j = 1:size(map,2)
%     y = map(j).rho_all;
%     p = map(j).sig_thresh;    
%     p = map(j).sig_thresh;    
%     for i = 1:10
%         temp = abs(y(:,:,i));
%         x(j,i) = 100*(sum(temp>p(i),'all')/sum(~isnan(temp(:))));        
%     end
% end
% % show that the first dimension spans the cortex
% figure; hold on;
% a = x(:,1);
% b = x(:,2:3);
% b(b<1)=[];
% b = b(:);
% c = x(:,4:10);
% c(b<1)=[];
% c = b(:);
% a = pairedBootstrap(a,@nanmean);
% b = pairedBootstrap(b,@nanmean);  
% c = pairedBootstrap(c,@nanmean); 
% CompareViolins(a',fp,'plotspread',0,'divfactor',1.,'plotaverage',1,'col',{fp.c_ff},'distWidth',0.75,'xpos',1);
% CompareViolins(b(:)',fp,'plotspread',0,'divfactor',1.,'plotaverage',1,'distWidth',0.75,'xpos',2);    
% CompareViolins(c(:)',fp,'plotspread',0,'divfactor',1.,'plotaverage',1,'distWidth',0.75,'xpos',3);  
% xlim([0.5 3.5])
% ylim([48 80])
% ylabel({'Network size','(% of dorsal cortex)'})
% xlabel('Dimension')
% set(gca,'xtick',[1,2,3],'ytick',[0:20:100])
% fp.FormatAxes(gca); box on; grid on; 
% fp.FigureSizing(gcf,[3 2 4 5],[2 10 10 10])
% 
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'Size','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\BehavioralRelationship',0); close all


end