function [ccaval,rrrval,rho,deltaD] = rrrCCAoverlap(data_rrr,data_cca,verbose)
%Camden MacDowell - timeless
%Compares the overlap in subspaces discovered from rrr and cca

% data_cca = load('Mouse 331 Recording 1CCA_muaflag0_motif3.mat');
% data_rrr = load('Mouse 331 Recording 1RRR_muaflag0_motif3.mat');

%get the subspaces
good_idx = isSubspace(data_rrr, 0);

%get subspaces in rrr
ccaval = NaN(size(data_rrr,2),2);    %percent of invalid in rrr that were also invalid in cca
rrrval = NaN(size(data_rrr,2),2);    %percent of invalid in rrr that were also invalid in cca
rho = NaN(1,size(data_rrr,2));   %correlation in the number of subspaces 
deltaD = NaN(1,size(data_rrr,2));  %mean difference in number of dimensions in rrr vs cca
for cur_fit = 1:size(data_rrr,2)    

    %parse areas
    a_cca = data_cca(cur_fit).paired_areas;
    a_rrr = data_rrr(cur_fit).paired_areas;

    %dimensions of rrr
    n_rrr = cellfun(@(x) ModelSelect([ mean(x); std(x)/sqrt(size(x,1)) ], 1:size(x,2)), data_rrr(cur_fit).cvl_rrr,'UniformOutput',1);

    %minimum of rrr local dimensions
    d = min(cat(1,data_rrr(cur_fit).qOpt{:}),cat(1,data_rrr(cur_fit).qOpt_target{:}));

    %remove self indices
    self_idx = diff(a_rrr,[],2)==0;
    a_rrr(self_idx,:)=[]; 
    good_idx{cur_fit}(self_idx)=[]; 
    n_rrr(self_idx)=[];
    d(self_idx)=[];

    %loop through each subspace 
    n_cca = NaN(size(a_rrr,1),1);
    for i = 1:size(a_rrr,1)
        %get the cca (1) that corresponds to the rrr directioned subspaces (2)
        [~,idx_cca,~] = intersect(a_cca,cat(1,a_rrr(i,:),fliplr(a_rrr(i,:))),'rows');    
        %dimensions of cca 
        n_cca(i) = numel(data_cca(cur_fit).r{idx_cca});
    end

    %set any violations to 'no subspace'
    n_cca(n_cca>=d)=0;

    %# where no subspace on cca and yes on rrr
    ccaval(cur_fit,1) = sum(good_idx{cur_fit}==1 & n_cca==0)/sum(n_cca==0);
    ccaval(cur_fit,2) = sum(n_cca==0); 

    %how many nonsubspaces also show none in CCA
    rrrval(cur_fit,1) = sum(good_idx{cur_fit}==0 & n_cca==0)/(sum(good_idx{cur_fit}==0));
    rrrval(cur_fit,2) = sum(good_idx{cur_fit}==0); 

    %overall correlation in dimensions (after removing the non)
    bad_idx = (good_idx{cur_fit}==0 | n_cca==0);

    n_cca(bad_idx)=[];
    n_rrr(bad_idx)=[];
    rho(cur_fit) = corr(n_cca,n_rrr,'type','Spearman');
    deltaD(cur_fit) = nanmean(n_rrr-n_cca);
    
    if verbose
        figure('units','normalized','position',[ 0.3542    0.5167    0.5698    0.3889]); hold on; 
        rng('default');
        subplot(1,2,1); hold on;
        plot(n_cca+rand(numel(n_cca),1)/4,n_rrr+rand(numel(n_rrr),1)/4,'marker','o','MarkerEdgeColor',[0.5 0.5 0.5],'linestyle','none','linewidth',1.5)
        plot(1:10,1:10,'linestyle','--');
%         legend('strong','medium','weak','unity');
        set(gca,'xlim',[1,max(cat(1,n_cca,n_rrr))],'ylim',[1,max(cat(1,n_cca,n_rrr))])
        axis square;
        title(sprintf('rho %0.2f',rho(cur_fit)),'fontweight','normal')
        xlabel('CCA dimensions')
        ylabel('RRR dimensions')
        subplot(1,2,2); hold on;
        histogram(n_rrr-n_cca,'BinWidth',1,'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5)
        plot([0 0],get(gca,'ylim'),'color','k','linewidth',2,'linestyle','--')
        title('# RRR Dim - # CCA dim','Fontweight','normal')
        axis square;
    end
end

end %function end
















