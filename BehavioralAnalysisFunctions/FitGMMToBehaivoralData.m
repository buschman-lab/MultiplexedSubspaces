function [clust_prob,clust_id] = FitGMMToBehaivoralData(features_downsampled,verbose)
if nargin <2; verbose = 0; end
rng('default');
Y = features_downsampled;
options = statset('Display','final','MaxIter',500);
gmm = fitgmdist(Y,2,'options',options,'CovarianceType','full','RegularizationValue',0,'SharedCovariance',false);
[clust_idx,nlogL,P,logpdf,d2] = cluster(gmm,Y);

%smooth the probabilties. e.g. get a running product
kernel = 26;
P_smooth = convn(P,ones(kernel,1)/kernel,'same');

%get the max probabilty for each timepoint
[clust_prob,clust_id] = max(log(P_smooth),[],2);
% 
% % %remove points with low probabilties
% bad_id = clust_prob<=log(0.50);
% clust_prob(bad_id)=[];
% clust_id(bad_id)=[];
% Y(bad_id,:)=[];

if verbose
    %histogram of the full distribution and each sub_distribution
    num_states = numel(unique(clust_id));
    col = getColorPalet(num_states);
    figure('position',[680   387   458   604]); hold on; 
    for i = 1:size(Y,2)
        figure('position',[680   387   458   604]); hold on;
    %    subplot(3,1,i); hold on;
       histogram(Y(:,i),'BinWidth',0.1,'FaceAlpha',0.2,'FaceColor','k','EdgeColor','none')
       for j = 1:numel(unique(clust_id))
           histogram(Y(clust_id==j,i),'BinWidth',0.1,'EdgeColor','none','FaceColor',col(j,:),'FaceAlpha',0.7)
       end  
       xlim([-4 4])
       ylabel('Bins');
       if i ==size(Y,2)       
           xlabel('Z-Score')       
       end
       legend('location','EastOutside')
    %    setFigureDefaults;
       grid off
    %    set(gca,'position',[3 3 3 4])
    end
end
end