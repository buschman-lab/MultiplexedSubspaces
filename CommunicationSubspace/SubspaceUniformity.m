function [xboot,x,stats,ent] = SubspaceUniformity(B,AArea,AreaLabel,data,targ_area,cutoff,shufseed)
%% Get the fraction of top beta weights contributed by each area across motifs and recordings
if nargin <2; targ_area='VIS'; end
if nargin <3; cutoff=0.25; end
% if nargin <4; shufseed=1; end
[~, area_all] = LoadVariable(data,'rrr_beta',targ_area,1); %load betas
area_all(strcmp(area_all,targ_area))=[];

%get top x% of beta weights
nneu = NaN(6,14,7);
eachneu = NaN(6,14,7);
totneu = NaN(6,14);
for cur_rec = 1:6
    for cur_m = 1:14 
        area = AArea{cur_rec,cur_m};
        b = B{cur_rec,cur_m};
        area_label = AreaLabel{cur_rec,cur_m};
        %use absolute value
        b = abs(b); 
        [b,idx] = sort(b,'descend');
        area = area(idx);
        if shufseed>1 %shuffle labels            
           rng(shufseed);
           areatemp = area(randperm(numel(area),numel(area)));
        else
            areatemp = area;
        end
        totneu(cur_rec,cur_m) = numel(b);
        y = arrayfun(@(n) nansum(areatemp==n),unique(areatemp),'UniformOutput',1);
        idx = ceil(cutoff*numel(b));
        temp = arrayfun(@(n) nansum(areatemp(1:idx)==n),unique(areatemp),'UniformOutput',1);
        %adjust for recs missing locations
        area_miss = ismember(area_all,area_label);
        if ~isempty(temp)
            nneu(cur_rec,cur_m,area_miss)=temp;
            eachneu(cur_rec,cur_m,area_miss)=y;
        end
    end
end

%gutcheck that this should be around 20%
% n = 100*squeeze(nansum(nneu,3))./totneu; %this 

%Plot the relative contribution to betas normalized by the number
%of neurons in that area
x = 100*(nneu./eachneu);

%convert back to zeros where appropriate since that NaN value we don't want
%to lose. 
x(eachneu==0)=0;
% x = squeeze(nanmean(x,1));

x = reshape(x,size(x,1)*size(x,2),size(x,3));

% bootstrap across recordings and motifs
[xboot,stats] = pairedBootstrap(x,@nanmean);
stats.xboot = xboot; 

%get entropy of each motif/rec
get_ent = @(p) -sum(p.*log2(p)); 
ent = arrayfun(@(n) get_ent(x(n,:)/100), 1:size(x,1),'UniformOutput',1);
ent = nanmean(ent);

% % % %Shuffle x area labels 1000 times get entropy using the same base data
% % ent_null = ones(1,1000); 
% % rng('default');
% % for i = 1:1000    
% %     idx = arrayfun(@(n) randperm(7,7), 1:size(x,1),'UniformOutput',0);
% %     idx = cat(1,idx{:});
% %     temp = arrayfun(@(n) get_ent(x(n,:)/100), 1:size(x,1),'UniformOutput',1);    
% %     ent_null(i) = nanmean(temp);
% % %     p = nanmean(x(idx))/100;
% % %     ent_null(i) = -sum(p.*log2(p));
% % end




end %function end


function [b,area,area_label] = loadVisBeta(data,cur_d,cur_rec,cur_m,targ_area)
[b, area_label] = LoadVariable(data(cur_rec),'rrr_beta',targ_area,cur_d); %load betas
area_label = area_label(strcmp(area_label,targ_area)==0);
area=LoadVariable(data(cur_rec),'beta_region',targ_area);
%get one rec and one motif
b = (squeeze(b(cur_m,:)));
area = squeeze(area(cur_m,:));
idx = isnan(area);
area(idx)=[];
b(idx)=[];

%reorder by decreasing median beta
[b,reord] = sort(b,'descend');
area = area(reord);
end
