function Figure2d(x,ci,area_all)
fp = fig_params_cortdynamics;

idx = 0:0.1:1; %sweep fractions
%plot the chart
figure; hold on; 
col = fp.c_area; col = col(strcmp(area_all,area_all{8})==0,:);
for i = 1:7    
    tempci = abs(flipud(squeeze(ci(i,:,:))-x(i,:)));
    shadedErrorBar(idx,x(i,:),tempci,'lineprops',{'color',[col(i,:) 0.75],'linewidth',2},'patchSaturation',0.15);
end
ylabel('Percentage of Population');
xlabel('Fraction of Beta weights');
fp.FormatAxes(gca); box on; grid on
fp.FigureSizing(gcf,[3 3 4 4],[10 10 10 10])
title('Visual dimesion 1');    

%plot the AUC in a small panel
%mean auc
xx = x/100;
auc = NaN(1,size(x,1));
aucci = NaN(2,size(x,1));

for i = 1:size(x,1)
    auc(i) = trapz(xx(i,:)/numel(xx(i,:)));
    temp = squeeze(ci(i,1,:))/100;
    aucci(1,i) = trapz(temp/numel(temp));
    temp = squeeze(ci(i,2,:))/100;
    aucci(2,i) = trapz(temp/numel(temp));
end

%mini plot        
a = area_all(strcmp(area_all,area_all{8})==0);
[auc,idx] = sort(auc,'descend');
a = a(idx);
col = fp.c_area; col = col(strcmp(area_all,area_all{8})==0,:);
col = col(idx,:);
aucci = aucci(:,idx);
figure; hold on; 
for i = 1:numel(a)
    errorbar(i,auc(i),auc(i)-aucci(1,i),aucci(2,i)-auc(i),'linestyle','none','marker','o',...
        'MarkerEdgeColor','none','MarkerFaceColor',col(i,:),'markersize',2,'linewidth',1.5,'color',col(i,:),'CapSize',4)            
end
xlim([0 8])
plot([0,8],[0.5 0.5],'linestyle','--','color',[0.25 0.25 0.25],'linewidth',1);
ylim([0.35,0.6])
fp.FormatAxes(gca); box on; grid on
fp.FigureSizing(gcf,[3 3 1.5 1.5],[10 10 12 10]) 

end