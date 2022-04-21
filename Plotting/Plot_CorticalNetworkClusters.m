function num_clust = Plot_CorticalNetworkClusters(data,type,verbose)
if nargin <3; verbose = 0; end %plotting of similarity matrices
fp = fig_params_cortdynamics;
%Across all recordings and areas
rng('default')
knn = 12; 
lvr = 1;
num_clust = [];    
for cur_rec = 1:6
    switch type
        case 'overall'           
            x = arrayfun(@(n) conditionDffMat(data{cur_rec}(n).rho_all),1:size(data{cur_rec},2),'UniformOutput',0);
            x = cat(1,x{:});
            tcorr_mat = abs(corr(x'));
            [cluster_idx, ~, ~] = PhenoCluster(tcorr_mat,'k',knn,'louvain_restarts',lvr,'verbose',0);           
            cluster_idx = cluster_idx(:,end);
            num_clust(cur_rec) = numel(unique(cluster_idx));
            %visualize cross correlation matrix
            if verbose
                fprintf('\n\tPlotting similarity matrix')
                Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx,'orderbydendrogram',1);
                title('Similarity Between Cortical Networks','fontsize',5,'fontweight','normal')
            end
        case 'wMotifxRegion'
            y = data{cur_rec};
            %loop through motifs
            m = unique(cat(1,y.cur_motif));
            for i = 1:numel(m)               
                idx = find(cat(1,y.cur_motif)==m(i)); %instances of motif
                x = arrayfun(@(n) conditionDffMat(y(n).rho_all(:,:,1:10)), idx,'UniformOutput',0);
                x = cat(1,x{:});
                tcorr_mat = abs(corr(x'));
                [cluster_idx, ~, ~] = PhenoCluster(tcorr_mat,'k',knn,'louvain_restarts',lvr,'verbose',0);           
                cluster_idx = cluster_idx(:,end);
                num_clust(cur_rec,i) = numel(unique(cluster_idx));
                %visualize cross correlation matrix     
                if verbose
                    fprintf('\n\tPlotting similarity matrix')
                    Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx,'orderbydendrogram',1);
                    title(sprintf('wMotifxRegion M | %d',i),'fontsize',12,'fontweight','normal')
                end
            end
        case 'wRegionxMotif'
            y = data{cur_rec};
            %loop through motifs
            m = unique(cat(1,y.cur_a));
            for i = 1:numel(m)               
                idx = find(cat(1,y.cur_a)==m(i)); %instances of motif
                x = arrayfun(@(n) conditionDffMat(y(n).rho_all), idx,'UniformOutput',0);
                x = cat(1,x{:});
                tcorr_mat = abs(corr(x'));
                [cluster_idx, ~, ~] = PhenoCluster(tcorr_mat,'k',knn,'louvain_restarts',lvr,'verbose',0);           
                cluster_idx = cluster_idx(:,end);
                num_clust(cur_rec,i) = numel(unique(cluster_idx));
                %visualize cross correlation matrix
                if verbose
                    fprintf('\n\tPlotting similarity matrix')
                    Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx,'orderbydendrogram',1);
                    title(sprintf('wRegionxMotif A | %d',cur_a),'fontsize',12,'fontweight','normal')
                end                                               
            end
            if cur_rec == 3
                temp = num_clust(cur_rec,:);
                temp(end-1:end)=[];
                num_clust(cur_rec,:)=NaN;
                num_clust(cur_rec,[1,2,3,5,7,8])=temp;
            elseif cur_rec==4
                temp = num_clust(cur_rec,:);
                temp(end)=[];
                num_clust(cur_rec,:)=NaN;
                num_clust(cur_rec,[1,2,3,4,5,7,8])=temp;
            end

        otherwise
            error('unknown clustering type');
    end
end


figure; hold on; 
histogram(num_clust(:),'BinWidth',1,'FaceAlpha',0.5,'FaceColor',[0.25 0.25 0.25],'EdgeAlpha',1,'EdgeColor',[0.25 0.25 0.25])
yval = get(gca,'ylim');
xlim([0 max(num_clust(:))+1])
plot([nanmean(num_clust(:)),nanmean(num_clust(:))],yval,'linewidth',2,'color','k','linestyle','-')
xlabel('rho');
ylabel('# of trial pairs');
fp.FormatAxes(gca); 
box on;
title({'Number of Clusters',sprintf('r=%0.2f',nanmean(num_clust(:)))},'FontWeight','normal')
fp.FigureSizing(gcf,[5 2 4 4],[2 10 12 12])

end

%to plot
%the number of clusters per motif across recordings
%the number of clusters per area across recordings
















