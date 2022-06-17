function PlotPowerLawDimensionality(data,cur_rec) 
fp = fig_params_cortdynamics;
%Camden 
ndim = 30;
%load data
[svca_mdl,~] = LoadVariable(data,'svca_across',[]);
svca_mdl = svca_mdl(:,:,:,1:ndim);
[rrr_mdl,area_label] = LoadVariable(data,'rel_performance',[],1);
rrr_mdl = rrr_mdl(:,:,:,1:ndim);
%plot an example of the power law decay in local and subspace dimensionality
plawplot(svca_mdl,rrr_mdl,area_label,cur_rec)

%then compare two histograms of the exponent (paired)
%get all exponents (for each recording)
local = getCoef(svca_mdl);
subspace = getCoef(rrr_mdl);

figure; hold on;
histogram(local(:),'binwidth',0.15,'Edgecolor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);
histogram(subspace(:),'binwidth',0.15,'Edgecolor',fp.c_ff,'FaceColor',fp.c_ff,'FaceAlpha',0.5);
xlabel('exponent')
ylabel({'models'});
[~,p] = ttest(local(:),subspace(:),'Tail','left');
title(sprintf('p=%0.3f',p),'FontWeight','normal');
fp.FormatAxes(gca); box on; grid off
xlim([0.25 3.]);
set(gca,'xtick',[0.5:1:2.5])
fp.FigureSizing(gcf,[3 2 4.5 4.5],[10 10 10 10])
col = fp.c_area;
Plot_CorrelateValuesBetweenRecordings(subspace',local','combo',fp,'right','xlabel',{'Subspace exponent'},...
    'ylabel','Local exponent','color_flag',1,'corrtype','kendall','col',col,'addjitter',3)
plot([0 max(cat(1,subspace(:),local(:)))],[0 max(cat(1,subspace(:),local(:)))],'linestyle',':','color','k')
xlim([0 2.5])
ylim([0 2.5])
title({'Subspace dimensionality decays','faster than local dim'},'fontweight','normal') 
fp.FigureSizing(gcf,[3 2 4.5 4.5],[10 10 10 10])


% plot by area
%reorder 
x = subspace./local;
[xboot,stats] = pairedBootstrap(x,@nanmean);
[~,idx] = sort(nanmean(x),'descend');
xboot = xboot(:,idx);   
% xboot = x(:,idx);

%get full stats accross them all
%     [~,statsall] = pairedBootstrap(x(:),@nanmean);

% plot as a violin showing each region
fp = fig_params_cortdynamics;
col = fp.c_area;
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
col = col(idx);
area_name = area_label(idx);

%flatten per area
figure; hold on;
vp = CompareViolins(xboot',fp,'col',col,'connectline',[],'plotspread',0,'divfactor',.25);
fp.FigureSizing(gcf,[3 2 6.5 4],[10 10 15 10]); 
set(gca,'XTickLabel',area_name,'XTickLabelRotation',45)
fp.FormatAxes(gca); box on; grid on
ylim([0.75 ceil(max(get(gca,'ylim')))])
plot(get(gca,'xlim'),[1 1],'linestyle','--','color','r','linewidth',2);
ylabel('Sub D / Local D');

%plot just showing the powerlaw dimensionality for local

x = local;
[xboot,stats] = pairedBootstrap(x,@nanmean);
[~,idx] = sort(nanmean(x),'descend');
xboot = xboot(:,idx);   
% xboot = x(:,idx);

%get full stats accross them all
%     [~,statsall] = pairedBootstrap(x(:),@nanmean);

% plot as a violin showing each region
fp = fig_params_cortdynamics;
col = fp.c_area;
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
col = col(idx);
area_name = area_label(idx);

%flatten per area
figure; hold on;
vp = CompareViolins(xboot',fp,'col',col,'connectline',[],'plotspread',0,'divfactor',.25);
fp.FigureSizing(gcf,[3 2 6.5 4],[10 10 15 10]); 
set(gca,'XTickLabel',area_name,'XTickLabelRotation',45)
fp.FormatAxes(gca); box on; grid on
% ylim([1.5 ceil(max(get(gca,'ylim')))])
% plot(get(gca,'xlim'),[1 1],'linestyle','--','color','k');
ylabel('Local D exponent');


%plot just showing the powerlaw inter-regional
x = subspace;
[xboot,stats] = pairedBootstrap(x,@nanmean);
[~,idx] = sort(nanmean(x),'descend');
xboot = xboot(:,idx); 
% plot as a violin showing each region
fp = fig_params_cortdynamics;
col = fp.c_area;
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
col = col(idx);
area_name = area_label(idx);

%flatten per area
figure; hold on;
vp = CompareViolins(xboot',fp,'col',col,'connectline',[],'plotspread',0,'divfactor',.25);
fp.FigureSizing(gcf,[3 2 6.5 4],[10 10 15 10]); 
set(gca,'XTickLabel',area_name,'XTickLabelRotation',45)
fp.FormatAxes(gca); box on; grid on
% ylim([1.5 ceil(max(get(gca,'ylim')))])
% plot(get(gca,'xlim'),[1 1],'linestyle','--','color','k');
ylabel('Interregional D exponent');


end %function end


function plawplot(svca_mdl,rrr_mdl,area_label,cur_rec)
rng('default')
fp = fig_params_cortdynamics;
warning off
%plot across all models for a region 
for i = 1:8
    x = squeeze(rrr_mdl(:,:,i,:));
    x = cat(3,x(:,:,1),diff(x,[],3));
%     x = squeeze(x(cur_rec,:,:)); %uncomment to just do one rec


    y = squeeze(svca_mdl(:,:,i,:));
    y = cat(3,y(:,:,1),diff(y,[],3));
%     y = squeeze(y(cur_rec,:,:)); %uncomment to just do one rec

    if size(x,3)>1
        x = reshape(x,size(x,1)*size(x,2),size(x,3));    
        y = reshape(y,size(y,1)*size(y,2),size(y,3));
    end       
    
    col = fp.c_area;
    figure; hold on;
    shadedErrorBar(1:size(x,2),nanmean(x),sem(x,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});
    shadedErrorBar(1:size(y,2),nanmean(y),sem(y,1),'lineprops',{'color',[0.1 0.1 0.1,0.50],'linewidth',2});
    set(gca,'xscale','log','yscale','log')
    xcoef = polyfit(log(1:size(x,2)),log(nanmean(x)),1); %slope, intercept
%     f = @(x) abs(xcoef(2))*(x^xcoef(1)); 
    f = @(x) (x^xcoef(1)); 
    fplot(f,[1,size(x,2)],'LineStyle','-','linewidth',2,'color',col(i,:),'linestyle','--');
    coef = polyfit(log(1:size(y,2)),log(nanmean(y)),1); %slope, intercept
%     f = @(x) abs(coef(2))*(x^coef(1)); 
    f = @(x) (x^coef(1)); 
    fplot(f,[1,size(y,2)],'LineStyle','-','linewidth',2,'color',[0.5 0.5 0.5],'linestyle','--');   

    xall = getAllCoef(x);
    yall = getAllCoef(y);
    cix = bootci(1000,@nanmean,xall);
    ciy = bootci(1000,@nanmean,yall);

    xlabel('dimension')
    ylabel({'variance explained'});
    title(sprintf('rec%d %s | 1/n^{%0.2f} CI: %0.2f %0.2f \n 1/n^{%0.2f} CI: %0.2f %0.2f ',cur_rec,area_label{i},nanmean(xall),cix(1),cix(2),nanmean(yall),ciy(1),ciy(2)),'fontsize',16,'fontweight','normal')
    fp.FormatAxes(gca); box on; 
    fp.FigureSizing(gcf,[3 2 3.5 3.5],[10 10 10 10])
end
warning on
end

function y = getCoef(xx)
y = NaN(size(xx,1)*size(xx,2),8);
warning off
for i = 1:8
    x = squeeze(xx(:,:,i,:));            
    x = cat(3,x(:,:,1),diff(x,[],3));
    
    x = reshape(x,size(x,1)*size(x,2),size(x,3));   
    
    y(:,i) = getAllCoef(x);
end


%%uncomment below for average by day
% y = NaN(6,8);
% for cur_rec = 1:6
%     for i = 1:8
%         x = squeeze(xx(:,:,i,:));
%         x = cat(3,x(:,:,1),diff(x,[],3));
%         x = squeeze(x(cur_rec,:,:));
% 
%         if size(x,3)>1
%             x = reshape(x,size(x,1)*size(x,2),size(x,3));    
%         end
% 
%         xcoef = polyfit(log(1:size(x,2)),log(nanmean(x)),1); %slope, intercept
%         y(cur_rec,i)=abs(real(xcoef(1)));
%     end
% end
end

function y = getAllCoef(x)
y = NaN(size(x,1),1);
for i = 1:size(x,1)
    if sum(isnan(x(i,:)))>1
        xtemp = x(i,1:find(isnan(x(i,:)),1,'first')-1); %some may not be all 30 dimensions
    else
        xtemp = x(i,:);
    end
    
    xcoef = polyfit(log(1:size(xtemp,2)),log(xtemp),1); %slope, intercept
    y(i) = abs(real(xcoef(1)));
    
    %deal with NaN filled recs missing certain areas
    if numel(xtemp)==0
        y(i)= NaN;
    end
end
end














