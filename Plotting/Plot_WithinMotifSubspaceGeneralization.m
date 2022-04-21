function Plot_WithinMotifSubspaceGeneralization(data)
%Camden MacDowell - timeless
%input data is the paired data
fp = fig_params_cortdynamics;

%% test an example motif
% get the coefficients for two places it projects to 
cur_motif = 5; 
r_self = cell(1,6);
r_null = cell(1,6);
r = cell(1,6);
r_max = cell(1,6);
for cur_rec = 1:6
    x = loadFunc(data,cur_rec,cur_motif,'VIS'); 
    y = loadFunc(data,cur_rec,cur_motif,'HIPP');
    [~,~,lambda] = loadCoef(data,cur_rec,cur_motif,{'VIS','PRE'}); %use the original ridge lambda
    [B,~] = loadCoef(data,cur_rec,cur_motif,{'VIS','PRE'});
    [BB,~] = loadCoef(data,cur_rec,cur_motif,{'VIS','HIPP'});
    % get correlation in betas
    rho = arrayfun(@(n) abs(corr(x*B(:,n),x*BB(:,n))),1:min(size(B,2),size(BB,2)),'UniformOutput',1); %projection area 1 to 2 using coef for 1
    rho_max = arrayfun(@(n) max(abs(corr(x*B(:,n),x*BB))),1:min(size(B,2),size(BB,2)),'UniformOutput',1); %projection area 1 to 2 using coef for 1
    rho_self = InternalXValidation(x,y,lambda);
    rng('default');
    rho_null = NaN(1,1000);
    for i = 1:1000        
        rho_null(i) = abs(corr(x*B(:,randperm(10,1)),x*BB(:,randperm(10,1))));
    end

    r_self{cur_rec} = fisherInverse(nanmean(fisherZ(rho_self)));
    r_null{cur_rec} = repmat(fisherInverse(nanmean(fisherZ(rho_null))),1,10);
    r_max{cur_rec} = rho_max(1:10);
    r{cur_rec} = rho(1:10);
end
r_self = cat(1,r_self{:});
r_null = cat(1,r_null{:});
r = cat(1,r{:});
r_max = cat(1,r_max{:});

%% plot the example
figure; hold on; 
arrayfun(@(n) plot(r(n,:),'color',[0.8 0.1 0.1],'LineWidth',0.25),1:size(r,1)); %plot all
plot(fisherInverse(nanmean(fisherZ(r))),'color',[0.8 0.1 0.1],'LineWidth',2) %plot average
arrayfun(@(n) plot(r_null(n,:),'color',[0.5 0.5 0.5],'LineWidth',0.25,'linestyle','-'), 1:size(r,1));
plot(fisherInverse(nanmean(fisherZ(r_null))),'color',[0.5 0.5 0.5],'LineWidth',2) %plot average
arrayfun(@(n) plot(r_self(n,:),'color',[0.1 0.1 0.8],'LineWidth',0.25,'LineStyle','-'), 1:size(r,1));
plot(fisherInverse(nanmean(fisherZ(r_self))),'color',[0.1 0.1 0.8],'LineWidth',2) %plot average
arrayfun(@(n) plot(r_max(n,:),'color',[0.1 0.8 0.1],'LineWidth',0.25,'LineStyle','-'), 1:size(r,1));
plot(fisherInverse(nanmean(fisherZ(r_max))),'color',[0.1 0.8 0.1],'LineWidth',2) %plot average
xlabel('Subspace dimension');
ylabel('Similarity of \beta');
title(sprintf('Subspaces span\nmultiple areas'),'Fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])


% shadedErrorBar(1:10,nanmean(r),sem(r,1),'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
% shadedErrorBar(1:10,nanmean(r_null),sem(r_null,1),'lineprops',{'color',[0.5 0.5 0.5 0.75],'linewidth',2, 'linestyle',':'});
% shadedErrorBar(1:10,nanmean(r_self),sem(r_self,1),'lineprops',{'color',[0.5 0.5 0.5 0.75],'linewidth',2});


end %function end




function [B,V,lambda] = loadCoef(data,cur_rec,cur_motif,area_name) 
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    paired_areas = data{cur_rec}(cur_motif).paired_areas;
    
    a = find(strcmp(area_label,area_name{1}));
    b = find(strcmp(area_label,area_name{2}));
    idx = find(ismember(paired_areas,[a,b],'rows')); 
    
    B = data{cur_rec}(cur_motif).rrr_B{idx};
    V = data{cur_rec}(cur_motif).rrr_V{idx};
    
    %get the originally used lambda
    x = area_val{strcmp(area_label,area_name{1})};
    if size(x,3)>1
        x = x-nanmean(x,3);
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';        
    end
    dMaxShrink = .5:.01:1;
    lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);
    cvl_ridge = data{cur_rec}(cur_motif).cvl_ridge{idx};
    [~,idx] = bestLambda(cvl_ridge);
    lambda = lambda(idx);
    
end


function x = loadFunc(data,cur_rec,cur_motif,area_name,flag)
    if nargin <5; flag = 1; end
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    x = area_val{strcmp(area_label,area_name)};

    %normalize to baseline
    x = normalizeToBaseline(x,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    
    if flag
        %subtract the psth
        x = x-nanmean(x,3);

        %concatentate across trials and pca
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    end
end
