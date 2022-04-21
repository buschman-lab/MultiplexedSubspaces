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
xlim([0.25 2.75]);
set(gca,'xtick',[0.5:1:2.5])
fp.FigureSizing(gcf,[3 2 3.25 3.25],[10 10 10 10])



end %function end


function plawplot(svca_mdl,rrr_mdl,area_label,cur_rec)
rng('default')
fp = fig_params_cortdynamics;
warning off
for i = 1:8
    x = squeeze(rrr_mdl(:,:,i,:));
    x = cat(3,x(:,:,1),diff(x,[],3));
    x = squeeze(x(cur_rec,:,:));


    y = squeeze(svca_mdl(:,:,i,:));
    y = cat(3,y(:,:,1),diff(y,[],3));
    y = squeeze(y(cur_rec,:,:));
    
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
    f = @(x) 1/(x^abs(xcoef(1))); 
    fplot(f,[1,size(x,2)],'LineStyle','-','linewidth',2,'color',col(i,:),'linestyle','--');
    coef = polyfit(log(1:size(y,2)),log(nanmean(y)),1); %slope, intercept
    f = @(y) 1/(y^abs(coef(1))); 
    fplot(f,[1,size(y,2)],'LineStyle','-','linewidth',2,'color',[0.5 0.5 0.5],'linestyle','--');   
    
    xlabel('dimension')
    ylabel({'variance explained'});
    title(sprintf('rec%d %s | 1/n^{%0.2f} | 1/n^{%0.2f}',cur_rec,area_label{i},abs(xcoef(1)),abs(coef(1))),'fontsize',16,'fontweight','normal')
    fp.FormatAxes(gca); box on; 
    fp.FigureSizing(gcf,[3 2 3.25 3.25],[10 10 10 10])
end
warning on
end

function y = getCoef(xx)
y = NaN(6,8);
for cur_rec = 1:6
    for i = 1:8
        x = squeeze(xx(:,:,i,:));
        x = cat(3,x(:,:,1),diff(x,[],3));
        x = squeeze(x(cur_rec,:,:));

        if size(x,3)>1
            x = reshape(x,size(x,1)*size(x,2),size(x,3));    
        end

        xcoef = polyfit(log(1:size(x,2)),log(nanmean(x)),1); %slope, intercept
        y(cur_rec,i)=abs(real(xcoef(1)));
    end
end
end















