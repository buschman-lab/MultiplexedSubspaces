function [x,area_label,area,b] = getBeta(data,cur_d,cur_rec,cur_m,targ_area,source_area)
    if nargin <6; source_area = []; end
    [~,area_all] = LoadVariable(data(1),'rrr_beta',targ_area,cur_d); %load betas
    [b, area_label] = LoadVariable(data(cur_rec),'rrr_beta',targ_area,cur_d); %load betas
    area_label = area_label(strcmp(area_label,targ_area)==0);
    area=LoadVariable(data(cur_rec),'beta_region',targ_area);
    %get one rec and one motif
    b = (squeeze(b(cur_m,:)));
    area = squeeze(area(cur_m,:));
    idx = isnan(area);
    area(idx)=[];
    b(idx)=[];    
    
    %% just get the beta weights for that region
    if ~isempty(source_area)
        x = b(area==find(strcmp(area_all,source_area)))';    
    else
        x = [];
    end
end

% % Plot and example: HIPP or RSP to VIS and SS boostrap confidence intervals and significance
% x = getBeta(data,1,2,1,'PRE','HIPP');
% y = getBeta(data,1,2,1,'VIS','HIPP');
% figure; hold on; 
% plot(x,y,'marker','.','color','k','markersize',fp.markersizebig*1.25,'LineStyle','none');
% AddLSline(x,y,x,'k');
% rho = corr(x,y);