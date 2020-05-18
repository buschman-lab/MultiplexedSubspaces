function [kval,fh] = FitPhenoK(X,varargin)
%Camden MacDowell - timeless. Randomly resamples (with replacement)
%n x n similarity matrix X to determine the best K value for phenographs K
%nearest neighbor search. 

opts.k_range = [2,4,6,8,10,15,20,30,50];
opts.louvain_restarts =5;
opts.verbose = 1;
opts.number_resamples = 25;
opts.genfigs = 1; 
opts = ParseOptionalInputs(opts,varargin);

%for reproducibility
rng('default');

%catch to only use k values at smaller than 1/4th the the size of the data
opts.k_range(opts.k_range>floor(size(X,1))/4)=[]; 

ovr_q = NaN(opts.number_resamples,numel(opts.k_range));
num_clust = NaN(opts.number_resamples,numel(opts.k_range));
for i = 1:numel(opts.k_range) 
    if opts.verbose; fprintf('\n Completed k value %d of %d', i, numel(opts.k_range)); end
    for j = 1:opts.number_resamples
        y = datasample(1:1:size(X,1),size(X,1),'Replace',true);
        [idx, ~, ovr_q(j,i)] = PhenoCluster(X(y,y),'k',opts.k_range(i),'louvain_restarts',opts.louvain_restarts,'Verbose',0);   
        if size(idx,2)>1
            idx = idx(:,end); %phenograph likes to over cluster so use the least clustered level if multiple
        end
        num_clust(j,i) = numel(unique(idx));
    end
end

%normalize
avg_q = (nanmean(ovr_q)-min(nanmean(ovr_q)))/(max(nanmean(ovr_q))-min(nanmean(ovr_q)));
var_q = (nanvar(ovr_q)-min(nanvar(ovr_q)))/(max(nanvar(ovr_q))-min(nanvar(ovr_q)));

%interpolate and find intersection (where crosses zero). could also do this by fitting a polynomial
xq = linspace(opts.k_range(1),opts.k_range(end),100); 
avg_q_int = interp1(opts.k_range,avg_q,xq,'linear');
var_q_int = interp1(opts.k_range,var_q,xq,'linear');

%approximate intersection
idx = find(avg_q_int-var_q_int<0,1,'first');
kval = round(xq(idx),0); %k must be an integer

%optionally make figures;
if opts.genfigs
   figure; hold on; 
   plot(opts.k_range,avg_q,'o-','color','r','linewidth',1.25); 
   plot(opts.k_range,var_q,'o-','color','b','linewidth',1.25); 
   plot(kval,var_q_int(idx),'kx','MarkerSize',15,'linewidth',2)
   legend('Modularity','Variability','K value')
   ylabel('AU (normalized 0-1)');
   xlabel('Number of Neighbors');
   setFigureDefaults;
   fh = gcf;
else
   fh = [];
end

