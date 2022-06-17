function CorticalNetworks(cur_rec,cur_motif,cur_a,reverseFlag)
%Camden MacDowell - timeless


if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CorticalNetworks\';
    FitPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace\';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CorticalNetworks/';
    FitPath = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CommunicationSubspace/';
end

fprintf('\n\tWorking on area %d rec %d motif %d',cur_a,cur_rec,cur_motif);

[rec_name,ImgPath,~,~,~] = LoadDataDirectories(cur_rec);

% Load motif fit
if reverseFlag ==1
    fn = [FitPath,rec_name,sprintf('regRRR_muaflag1_GROUPEDREVERSEmotif%d.mat',cur_motif)];
else
    fn = [FitPath,rec_name,sprintf('regRRR_muaflag1_GROUPEDmotif%d.mat',cur_motif)];
end
data = load(fn); 

% Load imaging data
load(ImgPath,'data_norm','nanpxs');
%temporally bin to match the binning on the ephys
dff = data_norm(:,1:2:end)+data_norm(:,2:2:end);
trig_dff = ParseByOnset(dff',[],data.motif_onset,data.win,cur_motif);

%remove the rare edge case where a motif begins at the start (no baseline)
trig_dff = RemoveEdgeTrials({trig_dff});
trig_dff = trig_dff{1};

fprintf('\n\tDone Loading data');

n_perm = 1000;
area_label = data.area_label;
alpha = 0.05; %threshold 
ndim = 10; 

rho_all = NaN(68,68,ndim);
sig_thresh = NaN(1,ndim);
rho_mc = NaN(ndim,n_perm);
for cur_d = 1:ndim    
    fprintf('\n\tWorking on dimension %d of %d',cur_d,ndim);
    B = data.rrr_B{cur_a}(:,cur_d);

    if reverseFlag
        x =data.area_val{strcmp(area_label,area_label{cur_a})==1};
    else
        x = cat(1,data.area_val{strcmp(area_label,area_label{cur_a})==0});
    end
    y = trig_dff;

    %normalize ephys to baseline
    x = normalizeToBaseline(x,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    y = y(:,3:end,:);

    %subtract the psth
    x = x-nanmean(x,3);
    y = y-nanmean(y,3);   

    %save for permutations
    xx = x; 
    
    %concatentate across trials and pca
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
    

    %organize B so that most neurons are positively contributing positively
    if sum(B>0)<sum(B<0) 
        B = B*-1;
    end
    
    %get projection
    xhat = x*B;
    
%     %To plot an example for one pixel
%     figure; hold on; 
%     plot(xhat(:),y(:,1102),'marker','.','color',[0.15 0.15 0.15],'markersize',3,'linestyle','none')
%     AddLSline(xhat(:),y(:,1102),xhat(:),[0.8 0 0]);
%     fp = fig_params_cortdynamics;
%     ylabel('Subspace Activity');
%     xlabel({'Deconvolved FR','of one pixel'});
%     fp.FormatAxes(gca); box on; grid on
%     fp.FigureSizing(gcf,[3 2 3.5 3.5],[2 10 10 10])
%     saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('ExampleSchematic'),pwd,0); close all

    
    %get cortical network    
    rho = arrayfun(@(n) corr(xhat,y(:,n),'type','pearson'), 1:size(y,2));
    n = size(rho,2);
    rho_all(:,:,cur_d) = (conditionDffMat(rho,nanpxs));      

    %trial-perm
    fprintf('running permutations');
    rng('default');
    rho_perm = NaN(n_perm,n);
    for cur_perm = 1:n_perm        
        x_temp = reshape(xx(:,:,randperm(size(xx,3),size(xx,3))),[size(xx,1),size(xx,2)*size(xx,3)])';
        x_temp = reshape(x_temp,[size(x_temp,1),size(x_temp,2)*size(x_temp,3)]);
        xhat = x_temp*B;
        rho_perm(cur_perm,:) = arrayfun(@(n) corr(xhat,y(:,n),'type','pearson'), 1:size(y,2));
    end

    %get the corrected threshold
    rho_mc(cur_d,:) = max(abs(rho_perm),[],2);
    sig_thresh(cur_d) = prctile(rho_mc(cur_d,:),100*(1-alpha));                    
end

%save off
if reverseFlag 
    fn = [savedir,rec_name,sprintf('REVERSEarea%dmotif%d.mat',cur_a,cur_motif)]; 
else
    fn = [savedir,rec_name,sprintf('area%dmotif%d.mat',cur_a,cur_motif)]; 
end
save(fn,'rho_mc','sig_thresh','alpha','rho_all','area_label','cur_motif','cur_a','reverseFlag');   
    
end %function 
 















%     %beta pseudo permutation
%     rng('default');
%     rho_perm = NaN(100,n);
%     for cur_perm = 1:100
%         shufidx = PermuteWithinArea(B,data,1,cur_a);
%         xhat = x*B(shufidx);
%         rho_perm(cur_perm,:) = arrayfun(@(n) corr(xhat,y(:,n),'type','pearson'), 1:size(y,2));
%     end


% %% cluster correction
% %get the cluster size
% [sz,pxlIdx,cSize] = clustSize(p,nanpxs,primaryAlpha);
% 
% %%
% %permutation
% rng('default');
% sz_perm = NaN(100,1);
% for cur_perm = 1:100
%     shufidx = PermuteWithinArea(B,data,1,cur_a);
%     xhat = x*B(shufidx);
%     [~,pp] = arrayfun(@(n) corr(xhat,y(:,n),'type','pearson'), 1:size(y,2));
%     sz_perm(cur_perm) = clustSize(pp,nanpxs,primaryAlpha);
% end
% %%
% prctile(sz_perm,95)
% 
% function [sz,pxlIdx,cSize] = clustSize(p,nanpxs,alpha)
%     BW = (conditionDffMat(p<alpha,nanpxs));      
%     BW(isnan(BW))=0;        
%     %get cluster correction 
%     CC = bwconncomp(BW,8);    
%     [cSize,idx] = sort(cellfun(@(x) numel(x),CC.PixelIdxList),'descend'); 
%     pxlIdx = CC.PixelIdxList(idx);
%     sz = cSize(1);
% end















