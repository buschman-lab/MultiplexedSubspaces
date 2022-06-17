function CommunicationSubspace_MotifTriggered(motif,cur_rec,muaflag,useCCA)
%Camden MacDowell - timeless
%Runs through a pipeline of CCA analyses for a given motif
%INPUTs
%win = [-5, 15] (default). The window around each motif spike to use.
%negative values are before onset. positive after. 
%EphysPath; the path of the ap_opts.mat file
%motif_fits; paths to the BasisMotifFits for a given mouse. 

%%Yes I know that this is horribly ugly code

if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    if useCCA<4
        savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace\';
    else
        savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace_Lagged\';
    end
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    if useCCA <4
        savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CommunicationSubspace/';
    else
        savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CommunicationSubspace_Lagged/';
    end
end
tic
win=[-2 11]; %hardcoded right now. 
%starting
fprintf('Working on rec %d motif %d',cur_rec,motif);
%% Gathering Data
[rec_name,~,~,EphysPath,motif_fits] = LoadDataDirectories(cur_rec);

%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',muaflag,'depth_type','probe'); 
% st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);
st_norm = st_mat;
st_norm = cellfun(@(x) x(1:2:end,:)+x(2:2:end,:),st_norm,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct
motif_onset = cellfun(@(x) floor(x/2), motif_onset,'UniformOutput',0);

if motif > numel(motif_onset) %get null periods that of window length    
    fprintf('\n\t Running on NULL')
    [~,trig_st] = ParseNullPeriods([],st_norm,motif_onset,win,10);
else %parse motif onsets
    fprintf('\n\t Running on Motifs')    
    [~,trig_st] = ParseByOnset([],st_norm,motif_onset,win,motif);
end

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'general');

%remove neurons that fire less than 0.5spike/sec on average across trials
[~, inactive_idx] = RemoveInactiveNeurons(area_val, 0.25/7.5);

%remove the rare edge case where a motif begins at the start (no baseline)
area_val = RemoveEdgeTrials(area_val);

%clean up areas %third input is the min # of spikes to keep area
[area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 

%% Main
if useCCA==1
    fprintf('USING CCA')
    % analyze all pairs of regions
    paired_areas = nchoosek(1:numel(area_label),2); 
    %preallocate variables
    a = cell(size(paired_areas,1),1);
    b = cell(size(paired_areas,1),1);
    U = cell(size(paired_areas,1),1);
    V = cell(size(paired_areas,1),1);
    r = cell(size(paired_areas,1),1);
    r_baseline = cell(size(paired_areas,1),1);
    pval = cell(size(paired_areas,1),1);
    pval_robust = cell(size(paired_areas,1),1);
    t = cell(size(paired_areas,1),1);
    best_idx = cell(size(paired_areas,1),1);
    aTheta_xv = cell(size(paired_areas,1),1);
    bTheta_xv = cell(size(paired_areas,1),1);
    xx_coef = cell(size(paired_areas,1),1);
    yy_coef = cell(size(paired_areas,1),1);
    a_pca = cell(size(paired_areas,1),1);
    b_pca = cell(size(paired_areas,1),1);
    for i = 1:size(paired_areas,1)
        fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,1));
        x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
        y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);

        [a{i},b{i},U{i},V{i},r{i},pval{i},~,aTheta_xv{i},bTheta_xv{i}] = significantCVs(x,y,0.01,0);
    end %subspace identification loop
    % save off data
    save([savedir,sprintf('%sCCA_muaflag%d_motif%d',rec_name,muaflag,motif)])
%     saveCurFigs(get(groot, 'Children'),{'-dpng'},[rec_name,'CCAstrength',sprintf('motif %d',motif)],savedir,0); close all
    fprintf('\ndone')
    toc
    

elseif useCCA==0 %rrr per area
    fprintf('USING RRR')
    paired_areas = combvec(1:numel(area_label),1:numel(area_label))';
    %preallocate variables
    qOpt = cell(size(paired_areas,1),1);
    qOpt_target = cell(size(paired_areas,1),1);
    cvl_fa = cell(size(paired_areas,1),1);
    rrr_B = cell(size(paired_areas,1),1);
    rrr_V = cell(size(paired_areas,1),1);
    cvl_rrr = cell(size(paired_areas,1),1);
    cvl_ridge = cell(size(paired_areas,1),1);
    cvl_fa_target = cell(size(paired_areas,1),1);
    rrr_B_xval = cell(size(paired_areas,1),1);
    for i = 1:size(paired_areas,1)
        fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,1));
        x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
        y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);

        [qOpt{i},qOpt_target{i},cvl_fa{i},rrr_B{i},cvl_rrr{i},cvl_ridge{i},cvl_fa_target{i},rrr_B_xval{i},rrr_V{i}] = RRR(x,y);
    end %subspace identification loop
    % save off data
    save([savedir,sprintf('%sregRRR_muaflag%d_motif%d',rec_name,muaflag,motif)])
    fprintf('\ndone')
    toc
    
elseif useCCA==2
    fprintf('USING CCA_full')
    % analyze all pairs of regions
    paired_areas = 1:numel(area_label);
    %preallocate variables
    a = cell(size(paired_areas,1),1);
    b = cell(size(paired_areas,1),1);
    U = cell(size(paired_areas,1),1);
    V = cell(size(paired_areas,1),1);
    r = cell(size(paired_areas,1),1);
    pval = cell(size(paired_areas,1),1);
    aTheta_xv = cell(size(paired_areas,1),1);
    bTheta_xv = cell(size(paired_areas,1),1);
    for i = 1:size(paired_areas,2)
        fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
        idx = strcmp(area_label,area_label{paired_areas(i)});
        x = cat(1,area_val{idx==0,:});
        y = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);

        [a{i},b{i},U{i},V{i},r{i},pval{i},~,aTheta_xv{i},bTheta_xv{i}] = significantCVs_full(x,y,0.05,0);
    end %subspace identification loop
    % save off data
    save([savedir,sprintf('%sCCA_muaflag%d_GROUPEDmotif%d',rec_name,muaflag,motif)])
    fprintf('\ndone')
    toc
    
elseif useCCA==3 %rrr across areas
    fprintf('USING RRR')
    paired_areas = 1:numel(area_label);
    %preallocate variables
    rrr_b = cell(size(paired_areas,1),1);
    rrr_V = cell(size(paired_areas,1),1);
    rrr_B = cell(size(paired_areas,1),1);
    cvl_rrr = cell(size(paired_areas,1),1);
    rrr_B_xval = cell(size(paired_areas,1),1);
    cvl_ridge = cell(size(paired_areas,1),1);
    grouping = cell(size(paired_areas,1),1);
    cvl_rrr_unique = cell(size(paired_areas,1),1);
    contribution = cell(size(paired_areas,1),1);    
    for i = 1:size(paired_areas,2)
        fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
        idx = strcmp(area_label,area_label{paired_areas(i)});
        x = cat(1,area_val{idx==0,:});
        y = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);

        %while here, go ahead and parse beta by area
        grp = arrayfun(@(n) ones(size(area_val{n},1),1)*n,find(idx==0),'UniformOutput',0);
        grp = cat(1,grp{:});
        grouping{i} = grp;
        [rrr_b{i},rrr_B{i},rrr_V{i},cvl_rrr{i},cvl_ridge{i},cvl_rrr_unique{i},contribution{i},rrr_B_xval{i}] = RRR_full(x,y,grp);
    end %subspace identification loop
    % save off data
    save([savedir,sprintf('%sregRRR_muaflag%d_GROUPEDmotif%d',rec_name,muaflag,motif)])
%     save([savedir,sprintf('%sregRRR_muaflag%d_GROUPmeanmotif%d',rec_name,muaflag,motif)])
%     save([savedir,sprintf('%sregRRR_muaflag%d_GROUPnonemotif%d',rec_name,muaflag,motif)])
    fprintf('\ndone')
    toc    
elseif useCCA==4 %rrr across areas REVERSED
    fprintf('USING RRR')
    paired_areas = 1:numel(area_label);
    %preallocate variables
    rrr_b = cell(size(paired_areas,1),1);
    rrr_V = cell(size(paired_areas,1),1);
    rrr_B = cell(size(paired_areas,1),1);
    rrr_B_xval = cell(size(paired_areas,1),1);
    cvl_rrr = cell(size(paired_areas,1),1);
    cvl_ridge = cell(size(paired_areas,1),1);
    grouping = cell(size(paired_areas,1),1);
    cvl_rrr_unique = cell(size(paired_areas,1),1);
    contribution = cell(size(paired_areas,1),1);    
    for i = 1:size(paired_areas,2)
        fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
        idx = strcmp(area_label,area_label{paired_areas(i)});
        y = cat(1,area_val{idx==0,:});
        x = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);

        %while here, go ahead and parse beta by area
        grp = arrayfun(@(n) ones(size(area_val{n},1),1)*n,find(idx==0),'UniformOutput',0);
        grp = cat(1,grp{:});
        grouping{i} = grp;
        [rrr_b{i},rrr_B{i},rrr_V{i},cvl_rrr{i},cvl_ridge{i},cvl_rrr_unique{i},contribution{i},rrr_B_xval{i}] = RRR_full(x,y,grp);
    end %subspace identification loop
    % save off data
    save([savedir,sprintf('%sregRRR_muaflag%d_GROUPEDREVERSEmotif%d',rec_name,muaflag,motif)])
%     save([savedir,sprintf('%sregRRR_muaflag%d_GROUPmeanREVERSEmotif%d',rec_name,muaflag,motif)])
    fprintf('\ndone')
    toc   
elseif useCCA==5 %rrr across areas Reverse with Lag
    fprintf('USING RRR_withlag')
    paired_areas = 1:numel(area_label);
    %preallocate variables
    rrr_b = cell(size(paired_areas,1),1);
    rrr_V = cell(size(paired_areas,1),1);
    rrr_B = cell(size(paired_areas,1),1);
    rrr_B_xval = cell(size(paired_areas,1),1);
    cvl_rrr = cell(size(paired_areas,1),1);
    cvl_ridge = cell(size(paired_areas,1),1);
    grouping = cell(size(paired_areas,1),1);
    cvl_rrr_unique = cell(size(paired_areas,1),1);
    contribution = cell(size(paired_areas,1),1);    
    for i = 1:size(paired_areas,2)
        fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
        idx = strcmp(area_label,area_label{paired_areas(i)});
        y = cat(1,area_val{idx==0,:});
        x = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end-1,:);
        y = y(:,4:end,:);

        %while here, go ahead and parse beta by area
        grp = arrayfun(@(n) ones(size(area_val{n},1),1)*n,find(idx==0),'UniformOutput',0);
        grp = cat(1,grp{:});
        grouping{i} = grp;
        [rrr_b{i},rrr_B{i},rrr_V{i},cvl_rrr{i},cvl_ridge{i},cvl_rrr_unique{i},contribution{i},rrr_B_xval{i}] = RRR_full(x,y,grp);
    end %subspace identification loop
    % save off data
    save([savedir,sprintf('%sregRRR_muaflag%d_GROUPEDREVERSEmotif%d_lag1',rec_name,muaflag,motif)])
    fprintf('\ndone')
    toc    
elseif useCCA==6 %rrr across areas with Lag
    fprintf('USING RRR_withlag')
    paired_areas = 1:numel(area_label);
    %preallocate variables
    rrr_b = cell(size(paired_areas,1),1);
    rrr_V = cell(size(paired_areas,1),1);
    rrr_B = cell(size(paired_areas,1),1);
    rrr_B_xval = cell(size(paired_areas,1),1);
    cvl_rrr = cell(size(paired_areas,1),1);
    cvl_ridge = cell(size(paired_areas,1),1);
    grouping = cell(size(paired_areas,1),1);
    cvl_rrr_unique = cell(size(paired_areas,1),1);
    contribution = cell(size(paired_areas,1),1);    
    for i = 1:size(paired_areas,2)
        fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
        idx = strcmp(area_label,area_label{paired_areas(i)});
        x = cat(1,area_val{idx==0,:});
        y = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end-1,:);
        y = y(:,4:end,:);

        %while here, go ahead and parse beta by area
        grp = arrayfun(@(n) ones(size(area_val{n},1),1)*n,find(idx==0),'UniformOutput',0);
        grp = cat(1,grp{:});
        grouping{i} = grp;
        [rrr_b{i},rrr_B{i},rrr_V{i},cvl_rrr{i},cvl_ridge{i},cvl_rrr_unique{i},contribution{i},rrr_B_xval{i}] = RRR_full(x,y,grp);
    end %subspace identification loop
    % save off data
    save([savedir,sprintf('%sregRRR_muaflag%d_GROUPEDmotif%d_lag1',rec_name,muaflag,motif)])
    fprintf('\ndone')
    toc        
else
    
end %method if/else
    
end %function end



