function DetectChaos(cur_rec,muaflag,str,stp)
%Camden MacDowell
%shell for 'chaos.m' from toker et al., 2020 A simple method for detecting chaos in nature

if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\chaos'));
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\chaos\';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/chaos'));
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/analysisplayground/chaos/';
end

tic
fprintf('\n\t detecting chaos');
[rec_name,~,~,EphysPath,~] = LoadDataDirectories(cur_rec);

%load ephys data and downsample
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',muaflag,'depth_type','probe'); 
st_norm = st_mat;
st_norm = cellfun(@(x) x(1:2:end,:)+x(2:2:end,:),st_norm,'UniformOutput',0);

% get neural activity per area
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(2,st_norm{:})',neu_area,'general');

%remove neurons that fire less than 0.5spike/sec on average across trials
% area_val = cellfun(@(x) x(:,1:5000), area_val,'UniformOutput',0);
%grooming: str = 1; stp=675;
%rest, big pupil: str=3600; stp=4275; 
%rest, tiny pupil; str=5850; stp=6525; 
area_val = cellfun(@(x) x(:,str:stp), area_val,'UniformOutput',0);
[area_val, inactive_idx] = RemoveInactiveNeurons(area_val, 0.5/7.5);

%clean up areas %third input is the min # of spikes to keep area
[area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 

%reduce dim with PCA and compute chaos
rng('default')
c = cell(1,numel(area_val));
for i = 1:numel(area_val)
%    fprintf('\n\tworking on area %d of %d',i,numel(area_val));
%    [~, score, ~,~,latent] = pca(area_val{i}') ;
%    idx = find(cumsum(latent)>75,1,'first'); %take the first 90% of explained variance
%    %compute chaos
   fprintf('\n\tworking on area %d of %d',i,numel(area_val));
   temp = cell(1,size(area_val{i},1));
   for j = 1:size(area_val{i},1)
      fprintf('\n\tworking on neu %d of %d',j,size(area_val{i},1));
      try
         temp{j} = chaos(area_val{i}(j,:)); 
      catch 
      end
   end      
   c{i} = temp; 
end

save([savedir,sprintf('%sChaos_muaflag%d_str%d_stp%d',rec_name,muaflag,str,stp)])
fprintf('\ndone detecting chaos')
toc

end %function end 