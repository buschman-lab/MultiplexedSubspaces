%relationship to widefield

%get the Yhat for Dim 1 for each point in time and each trial and
%get the trial to trial variance in the top network. Then, per pixel,
%correlation the Y and Pixel value and get a resulting correlation map. 

%USE THE REVERSE | so it will be the subspace of vis to everywhere else
%could actually the opposite too and reinforce this reciprocity thing

%load example rec and data
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace\Mouse 331 Recording 1regRRR_muaflag1_GROUPEDmotif6.mat');

%% load the imaging data
[~,ImgPath,~,~,~] = LoadDataDirectories(1);
load(ImgPath,'data_norm','nanpxs');
%%
trig_dff = ParseByOnset(data_norm',[],motif_onset,win,motif);

%get y value for each timepoint
cur_a = 7;
figure; hold on; 
t=tiledlayout(3,3,'Padding','normal','TileSpacing','normal');
for cur_d = 1:9
    nexttile
    B = rrr_B{cur_a}(:,cur_d);
    V = rrr_V{cur_a}(:,cur_d);

    x = cat(1,area_val{strcmp(area_label,area_label{cur_a})==0});
    y = trig_dff;

    %normalize to baseline
    x = normalizeToBaseline(x,[1:2],'mean');
    y = normalizeToBaseline(y,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    y = y(:,3:end,:);

    %subtract the psth
    x = x-nanmean(x,3);
    y = y-nanmean(y,3);

    %concatentate across trials and pca
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

    xhat = x*B;

    %correlate the two 
    [rho,primaryP] = arrayfun(@(n) corr(xhat,y(:,n),'type','Spearman'), 1:size(y,2));
    n = size(rho,2);
    %threshold at p=0.05
    rho(primaryP>=0.05)=NaN;
    rho = (conditionDffMat(rho,nanpxs));    
    
    BW = primaryP<0.001;     
    BW = (conditionDffMat(primaryP<0.001,nanpxs));      
    BW(isnan(BW))=0;
    %get cluster correction 
    CC = bwconncomp(BW,4);
    [clustsize,idx] = max(cellfun(@(x) numel(x),CC.PixelIdxList));
    BW(CC.PixelIdxList{idx})=2;

    %cluster correction
    rng('default');
    clust_perm = NaN(100,1);
    for cur_perm = 1:100
       shufidx = PermuteWithinArea(B,data,1,cur_a);
       xhat = x*B(shufidx);
       [~,primaryP] = arrayfun(@(n) corr(xhat,y(:,n),'type','Spearman'), 1:size(y,2));
       BW = (conditionDffMat(primaryP<0.01,nanpxs));      
       BW(isnan(BW))=0;
       %get cluster correction 
       CC = bwconncomp(BW,4);
       [clust_perm(cur_perm),idx] = max(cellfun(@(x) numel(x),CC.PixelIdxList));     
    end
    
    %Significance threshold | correcting for multiple comparisons
    st = prctile(max(abs(rho_perm),[],2),95); %max correlation value for each permutation   

    imagesc(rho); colorbar; colormap magma 
    set(gca,'ydir','reverse')
    title(sprintf('%s IN | dimenion %d',area_label{cur_a},cur_d),'fontweight','normal');

end

%save off the thresholded and unthresholded figures

% set(gcf,'position',[450.0000  250.0000  832.6000  606.2000]);


