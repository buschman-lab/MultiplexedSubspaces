function CommunicationSubspace_MotifTriggered(motif,cur_rec)
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
tic
win=[-5 15]; %hardcoded write now. 
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

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%parse motif onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,win,motif);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%clean up areas %third input is the min # of spikes to keep area
[area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 

%% Main
%analyze all pairs of regions
paired_areas = nchoosek(1:numel(area_label),2); 
%preallocate variables
pca_theta = NaN(size(paired_areas));
a = cell(size(paired_areas,1),1);
b = cell(size(paired_areas,1),1);
U = cell(size(paired_areas,1),1);
V = cell(size(paired_areas,1),1);
r = cell(size(paired_areas,1),1);
pval = cell(size(paired_areas,1),1);
t = cell(size(paired_areas,1),1);
best_idx = cell(size(paired_areas,1),1);
aTheta_xv = cell(size(paired_areas,1),1);
bTheta_xv = cell(size(paired_areas,1),1);
for i = 1:size(paired_areas,1)
    fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,1));
    x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
    y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};
    
    %Identify any significant CV(subspaces) between each population
    [a{i},b{i},U{i},V{i},r{i},pval{i},t{i},best_idx{i},aTheta_xv{i},bTheta_xv{i}] = significantCVs_perTimePoint(x,y,0.05,0);
end %subspace identification loop


%to do: add the pca loop so done on the same time relationships as the CVs
% %pca loop
% for i = 1:size(paired_areas,1)
%     if ~isempty(r{i})
%         x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
%         y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};
%         %Get local activity pattern: the neural coefficients from PCA
%         [xw,yw] = localActivity(x,y,3);
%         
%         %Angle between local and shared variance 
%         pca_theta(i,1) = AngleBetweenWeights(xw,a{i},'none'); 
%         pca_theta(i,2) = AngleBetweenWeights(yw,b{i},'none');  
% 
%         %[pending] | percent variance captured by shared variance
%     end
% end %pca loop


%% save off data
save([savedir,sprintf('%s_newversion_motif%d',rec_name,motif)])
fprintf('\ndone')
toc

end 



