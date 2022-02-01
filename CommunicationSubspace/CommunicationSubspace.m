function CommunicationSubspace(motif,cur_rec)
%Camden MacDowell - timeless
%Runs through a pipeline of CCA analyses for a given motif
%INPUTs
%win = [-5, 15] (default). The window around each motif spike to use.
%negative values are before onset. positive after. 
%EphysPath; the path of the ap_opts.mat file
%motif_fits; paths to the BasisMotifFits for a given mouse. 

if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA\';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/analysisplayground/CCA/';
end

win=[0 15]; %hardcoded write now. 
%starting
fprintf('Working on motif %d',motif);

%% Gathering Data
[rec_name,~,~,EphysPath,motif_fits] = LoadDataDirectories(cur_rec);

%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',0,'depth_type','probe'); 
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(2,st_norm{:})',neu_area,'parent');

%clean up areas %third input is the min # of spikes to keep area
[area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 

%parse the weighings per motif
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%% Main 
%CCA on all pairs of regions
paired_areas = nchoosek(1:numel(area_label),2); 
%preallocate variables
a = cell(size(paired_areas,1),1);
b = cell(size(paired_areas,1),1);
U = cell(size(paired_areas,1),1);
V = cell(size(paired_areas,1),1);
r = cell(size(paired_areas,1),1);
r_norm = cell(size(paired_areas,1),1);
motif_cvs = cell(size(paired_areas,1),1);
for i = 1:size(paired_areas,1)
    fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,1));
    x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
    y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};    
    %Identify any significant CV(subspaces) between each population    
    [a{i},b{i},U{i},V{i},r{i},r_norm{i}] = significantCVs_entireRec(x,y,0.001);     
    
    %get the response of each subspace to motif onsets
    [~,trig_cv] = ParseByOnset([],{U{i},V{i}},motif_onset,win,motif);
    
    %get the motif that significant per trial
    motif_cvs{i} = getSigCV(trig_cv,0.001);
    
end %subspace identification loop

% save off data
save([savedir,sprintf('%s_fullrec_motif%d',rec_name,motif)])
fprintf('\ndone')

end 

function idx = getSigCV(trig_cv,alpha)
    n_perm=1000;
    %determine which cvs are significanctly expressed during this motif
    x = trig_cv{1}; 
    y = trig_cv{2};

    %average corr per trial
    r = NaN(size(x,3),size(x,1));
    for i = 1:size(x,3)
        temp_x = (x(:,:,i)-nanmean(x(:,:,i),2))';
        temp_y = (y(:,:,i)-nanmean(y(:,:,i),2))';
        r(i,:) = arrayfun(@(n) corr(temp_x(:,n),temp_y(:,n)),1:size(temp_x,2));
    end
    r = r.^2;

    %repeat 1000 times with trial-wise permuation to get 1000 averages
    r_perm = NaN(size(x,3),n_perm);
    rng('default');
    for j = 1:n_perm
        y_perm = y(:,:,randperm(size(y,3),size(y,3)));
        for i = 1:size(x,3)
            temp_x = (x(:,:,i)-nanmean(x(:,:,i),2))';
            temp_y = (y_perm(:,:,i)-nanmean(y_perm(:,:,i),2))';
            r_perm(i,j)  = corr(temp_x(:,1),temp_y(:,1));
        end
    end
    r_perm = r_perm.^2;

    %get pvalue for each
    r_perm_avg = nanmean(r_perm);
    r_avg = nanmean(r);
    pval = arrayfun(@(n) sum([r_perm_avg,r_avg(n)]>=r_avg(n))/(numel([r_perm_avg,r_avg(n)])),1:numel(r_avg));
    idx = find(pval<alpha);
end

