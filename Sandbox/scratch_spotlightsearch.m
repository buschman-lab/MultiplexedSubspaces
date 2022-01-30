%Spotlight search 

%Tonight to do: 
%get a list of the significant ones per brain region 
%run for all animals/recordings

%% SINGLE UNIT DECODING
sig_area = cell(1,6);
beta_all = cell(1,6);
sig_all = cell(1,6);
for cur_rec = 1:6
    cur_rec
    if cur_rec ==1
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==2
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_12_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==3
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==4
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_08_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_08_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==5
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==6
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_10_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    end
    motifs = [9,13];
    
    
    [auc_val,area_label,auc_null] = SingleUnitDecoder(EphysPath,motif_fits,motifs,[1,15], 'peak');


    rng('default');
    b = cellfun(@(x) nanmean(x), auc_val,'UniformOutput',1);
    b_null = cellfun(@(x) nanmean(x), auc_null,'UniformOutput',1);
    pval = NaN(1,numel(b)); 
    %both on sided tests
    for i = 1:numel(b)
        if b(i)>=0
            pval(i) = sum([b_null(i,:),b(i)]>=b(i))/numel([b_null(i,:),b(i)]);
        else
            pval(i) = sum([b_null(i,:),b(i)]<=b(i))/numel([b_null(i,:),b(i)]);
        end
    end

    %list of significant brain regions
    sig_area{cur_rec} = area_label(pval<0.05);
    beta_all{cur_rec} = b;
    sig_all{cur_rec} = pval;
    
    
    
end

%%

% create plots of the average decoding accuracy per brain region. Grey if
% significant with mulitple comparison correction. and colored for
% otherwise

%get all areas labels
area_label = cell(1,4);
for cur_rec = 1:4
    cur_rec
    if cur_rec ==1
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==2
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_12_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==3
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==4
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_08_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_08_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==5
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==6
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_10_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    end
    %load ephys data
    [st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe');
    st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);

    %get the anatomical locations
    neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
    neu_area = cat(2,neu_area{:});
    area_label{cur_rec} = unique(cat(1,neu_area.parent_label));
end

%unique areas
a = cat(1,area_label{:});
uniqa = unique(a);
b = cat(2,sig_all{:})';
c = cat(1,beta_all{:});


%remove the one neuron per recording from cc and fiber tracks
uniqa(strcmp(uniqa,'cc'))=[]; 
uniqa(strcmp(uniqa,'fiber tracts'))=[];
uniqa(strcmp(uniqa,'na'))=[];

col = getColorPalet(numel(uniqa));
figure; hold on; 
COUNT = 1;
rng('default');
m_list = {'x','^','*','d','s','h'};
m = arrayfun(@(n) repmat(m_list{n},numel(area_label{n}),1),1:numel(m_list),'uniformoutput',0);
m = cat(1,m{:});
%bar plot with dot plot on top
xtick_label = {};
for i = 1:numel(uniqa)
   idx = find(strcmp(a,uniqa(i))==1);
   if numel(idx)>=4
       bar(COUNT,nanmean(c(idx)),'facecolor',col(i,:),'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
       for j = 1:numel(idx)
           if b(idx(j))<0.05/20
               plot(COUNT-0.05+rand(1)/10,c(idx(j)),'color',[0.9 0.1 0.1],'marker',m(idx(j)),'markersize',6)
           else
               plot(COUNT-0.05+rand(1)/10,c(idx(j)),'color',[0.5 0.5 0.5],'marker',m(idx(j)),'markersize',6)
           end
       end
       xtick_label(COUNT) = uniqa(i);
       COUNT = COUNT+1;
   end
end
line(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle',':','linewidth',1)
ylim([0.45,0.6])
set(gca,'Xtick',1:numel(xtick_label),'XTickLabel',xtick_label)
title('Marker legend: M1=x/^ M2=*/d M3=s/h','fontweight','normal') 
ylabel('decoding AUC');


%% Population level

sig_area = cell(1,6);
beta_all = cell(1,6);
sig_all = cell(1,6);
stats = cell(1,6);
stats_null = cell(1,6);
for cur_rec = 1:6
    cur_rec
    if cur_rec ==1
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==2
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_12_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==3
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==4
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_08_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_08_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==5
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==6
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_10_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    end
    motifs = [9,13];
    
%     population level
    [beta,area_label{cur_rec},beta_null,stats{cur_rec},stats_null{cur_rec}] = PopulationDecoder(EphysPath,motif_fits,motifs,[1,15], 'peak');


    rng('default');
    b = cellfun(@(x) nanmean(x), beta,'UniformOutput',1);
    b_null = cellfun(@(x) nanmean(x), beta_null,'UniformOutput',1);
    pval = NaN(1,numel(b)); 
    %both on sided tests
    for i = 1:numel(b)
        if b(i)>=0
            pval(i) = sum([b_null(i,:),b(i)]>=b(i))/numel([b_null(i,:),b(i)]);
        else
            pval(i) = sum([b_null(i,:),b(i)]<=b(i))/numel([b_null(i,:),b(i)]);
        end
    end

    %list of significant brain regions
    sig_area{cur_rec} = area_label{cur_rec}(pval<0.05/20);
    beta_all{cur_rec} = b;
    sig_all{cur_rec} = pval;    
end

%unique areas
a = cat(1,area_label{:});
uniqa = unique(a);
b = cat(2,sig_all{:})';
c = cat(1,beta_all{:});


%remove the one neuron per recording from cc and fiber tracks
uniqa(strcmp(uniqa,'cc'))=[]; 
uniqa(strcmp(uniqa,'fiber tracts'))=[];
uniqa(strcmp(uniqa,'na'))=[];

m_list = {'x','^','*','d','s','h'};
m = arrayfun(@(n) repmat(m_list{n},numel(area_label{n}),1),1:numel(m_list),'uniformoutput',0);
m = cat(1,m{:});

col = getColorPalet(numel(uniqa));
figure; hold on; 
COUNT = 1;
rng('default');
%bar plot with dot plot on top
xtick_label = {};
for i = 1:numel(uniqa)
   idx = find(strcmp(a,uniqa(i))==1);
   if numel(idx)>=4
       bar(COUNT,nanmean(c(idx)),'facecolor',col(i,:),'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
       for j = 1:numel(idx)
           if b(idx(j))<0.05
               plot(COUNT-0.05+rand(1)/10,c(idx(j)),'color',[0.9 0.1 0.1],'marker',m(idx(j)),'markersize',6)
           else
               plot(COUNT-0.05+rand(1)/10,c(idx(j)),'color',[0.5 0.5 0.5],'marker',m(idx(j)),'markersize',6)
           end
       end
       xtick_label(COUNT) = uniqa(i);
       COUNT = COUNT+1;
   end
end
% ylim([0.45,0.6])
set(gca,'Xtick',1:numel(xtick_label),'XTickLabel',xtick_label)
title('Marker legend: M1=x/^ M2=*/d M3=s/h','fontweight','normal') 
ylabel('decoding AUC');


%% generate a quick snapshot of decoding ability across all motifs

%just get one directional motifs
motifs = []; COUNT=1;
for i = 1:15
    for j = 1:15
        if j>i
            motifs(COUNT,:) = [i,j];
            COUNT = COUNT+1;
        end
    end
end
motifs(sum([motifs(:,1)==2,motifs(:,2)==2],2)==1)=[]; %remove the noise motif

beta = cell(1,6);
area_label = cell(1,6);
auc_val = cell(1,6);
acc_val = cell(1,6);
for cur_rec = 1:6
    cur_rec
    if cur_rec ==1
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==2
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_12_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec ==3
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==4
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_08_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_08_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==5
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    elseif cur_rec==6
        EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_10_2021_RestingState_g1\ap_opts.mat';
        [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    end
    
    %population level
    [beta{cur_rec},area_label{cur_rec},auc_val{cur_rec},acc_val{cur_rec}] = PopulationDecoder(EphysPath,motif_fits,motifs,[1,15], 'peak',0);    
end

%get the average decodability between motifs across animals
all_auc = [auc_val{:}];
temp = NaN(15,15);
for i = 1:size(motifs)    
   temp(motifs(i,1),motifs(i,2)) = nanmean(all_auc(i,:));
end

close all; figure; hold on;
imagesc(temp,[0.5 1]); colormap magma
set(gca,'ydir','reverse','xdir','reverse','ylim',[0.5,15.5],'xlim',[0.5,15.5]);
for i = 1:size(motifs)    
   if sum(all_auc(i,:)>0.8)==6
%        rectangle('position',[motifs(i,2)-0.5,motifs(i,1)-0.5,1,1],'EdgeColor','g','linewidth',0.25);
       p1 = plot(motifs(i,2),motifs(i,1),'marker','d','color','g','linestyle','none');
   elseif sum(all_auc(i,:)>0.7)==6
%        rectangle('position',[motifs(i,2)-0.5,motifs(i,1)-0.5,1,1],'EdgeColor','c','linewidth',0.25);
       p2 = plot(motifs(i,2),motifs(i,1),'marker','o','color','c','linestyle','none');
%    elseif sum(all_auc(i,:)>0.6)==6
%        rectangle('position',[motifs(i,2)-0.5,motifs(i,1)-0.5,1,1],'EdgeColor','b','linewidth',0.25);
%        plot(motifs(i,2),motifs(i,1),'marker','o','color','b');
   elseif sum(all_auc(i,:)<0.5)>=1
       p3 = plot(motifs(i,2),motifs(i,1),'marker','x','color','w','linestyle','none');
%        rectangle('position',[motifs(i,2)-0.5,motifs(i,1)-0.5,1,1],'EdgeColor','w','linewidth',0.25);    
   else
   end
end

set(gca,'xtick',1:15,'ytick',1:15)

title('Decoding Accuracy | using all units | avg across recs/animals','Fontweight','normal');
xlabel('motif'); ylabel('motif');
c = colorbar;
ylabel(c,'rho');
legend([p1,p2,p3],{'= all >0.8','= all>0.7','= some<0.5'},'Location','southeast','Color','k','TextColor','w','Box','off');


%for comparison, get the spatial similarity of the motifs
W_basis = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\MotifDiscovery\Mouse_basis_motifs.mat','W_basis_Center');
W_basis = W_basis.W_basis_Center;
W_basis = nanmean(W_basis,3);
rho = triu(corr(W_basis),1);
figure; hold on; colormap magma
imagesc(rho,[0.5 1]); colorbar
set(gca,'ydir','reverse','xdir','reverse','ylim',[0.5,15.5],'xlim',[0.5,15.5]);
for i = 1:size(motifs)    
   if sum(all_auc(i,:)>0.8)==6
       p1 = plot(motifs(i,2),motifs(i,1),'marker','d','color','g','linestyle','none');
   elseif sum(all_auc(i,:)>0.7)==6
       p2 = plot(motifs(i,2),motifs(i,1),'marker','o','color','c','linestyle','none');
   elseif sum(all_auc(i,:)<0.5)>=1
       p3 = plot(motifs(i,2),motifs(i,1),'marker','x','color','b','linestyle','none');
   else
   end
end

title('Spatial Similarity Between Motifs','Fontweight','normal');
xlabel('motif'); ylabel('motif');
c = colorbar;
ylabel(c,'rho');
legend([p1,p2,p3],{'= all >0.8','= all>0.7','= some<0.5'},'Location','southeast','Color','k','TextColor','w','Box','off');
























