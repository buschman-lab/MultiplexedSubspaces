function [auc,pval,rr] = SubspaceROC_single(data,area_name,motifs,verbose)
%data is across recordings
%motifs is the two motifs to compare
%Camden MacDowell - timeless
if nargin <4; verbose =0;end
assert(numel(motifs)==2,'too many motifs')
ndim = 20; %some regions have fewer than 30; 
rr = NaN(numel(data),numel(motifs),ndim);
for i = 1:numel(data)
    count=1;
    for j = motifs
        idx = strcmp(data{i}(j).area_label,area_name);
        if sum(idx)>0
            %normalize to the max explained variance (full model)
            full_mdl = cellfun(@(x) nanmean(bestLambda(x,1)), data{i}(j).cvl_ridge(idx),'UniformOutput',1);
            %get the mean cross validated performance        
            temp = [0,(1-nanmean(data{i}(j).cvl_rrr{idx}(:,1:ndim)))/full_mdl,1];           
            rr(i,count,1:size(temp,2)) = temp;
        end
        count= count+1;
    end
end
%remove empty recs
bad_idx = nanvar(squeeze(rr(:,1,:)),[],2)<eps;
rr(bad_idx==1,:,:)=[];
a = squeeze(rr(:,1,:));
b = squeeze(rr(:,2,:));


%get the AUCs for each recording
auc = arrayfun(@(n) trapz(a(n,:),b(n,:)), 1:size(a,1));
if nanmean(auc) <0.5
    [~,pval] = ttest(auc,0.5,'tail','left');
else
    [~,pval] = ttest(auc,0.5,'tail','right');
end

if verbose
%plot across recordings
figure; hold on; 
plot(0:1,0:1,'linewidth',1,'color','r','linestyle',':');
arrayfun(@(n) plot(a',b','linewidth',0.25,'color',[0.75 0.75 0.75]),1:10)
plot(nanmean(a),nanmean(b),'linestyle','-','marker','none','linewidth',2,'color','k'); hold on; 
xlabel(sprintf('Motif %d',motifs(1)));
ylabel(sprintf('Motif %d',motifs(2)));
axis square
title(sprintf('auc=%0.2f, p=%0.2f',nanmean(auc),pval),'fontweight','bold')
end

end %function









