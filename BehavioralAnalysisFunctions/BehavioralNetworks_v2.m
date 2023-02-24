function BehavioralNetworks_v2(cur_rec,cur_motif,cur_a)
%Camden MacDowell - timeless
tic
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BehavioralNetworks\';
    FitPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace\';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/BehavioralNetworks/';
    FitPath = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CommunicationSubspace/';
end

fprintf('\n\tWorking on area %d rec %d motif %d',cur_a,cur_rec,cur_motif);

%load all the data and project along 
[rec_name,~,~,~,~] = LoadDataDirectories(cur_rec);

% Load motif fit
fn = [FitPath,rec_name,sprintf('regRRR_muaflag1_GROUPEDmotif%d.mat',cur_motif)];
data = load(fn); 

%get raw ephys data
[area_val, area_label] = ParseByArea(cat(2,data.st_norm{:})',data.neu_area,'general');
%clean up areas %third input is the min # of spikes to keep area
[area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 

% Load imaging data
[me, ~,~,~,~,~] = AnalyzeBehavior(cur_rec,0);
% y = cat(1,nanmean(abs(me(2:end,:))),abs(me(2:end,:))); 
y = nanmean(abs(me(2:end,:)));

fprintf('\n\tDone Loading data');

n_perm = 1000;
area_label = data.area_label;
alpha = 0.05; %threshold 
ndim = 10; 

rho_all = NaN(size(y,1),ndim);
lag_all = NaN(size(y,1),ndim);
rho_perm_all = NaN(size(y,1),ndim);
for cur_d = 1:ndim    
    fprintf('\n\tWorking on dimension %d of %d',cur_d,ndim);
    B = data.rrr_B{cur_a}(:,cur_d);
    
    %get the area of interest
    x = cat(1,area_val{strcmp(area_label,area_label{cur_a})==0});
    
    %get projection
    xhat = x'*B;
    
%     %visualize
%     figure; hold on; plot(xhat(1:10000),'linewidth',2,'color','k');
%     yyaxis right
%     plot(y(:,1:10000)');
%     figure; hold on; 
%     [temp,lag] = xcorr(xhat-nanmean(xhat),y-nanmean(y),15,'normalized');
%     plot(lag,temp);

    %correlate projection to behavior (doesn't have to be zero lag - because we are just interested in any relationship)
    [rho,lag] = getRho(xhat,y);
    
    %save off
    rho_all(:,cur_d)=rho;
    lag_all(:,cur_d)=lag;  
        
    %shuffle test (take the top value of this across all
    %dimensions and recordings - the 95% of that is the FDR correct sig)
    rng('default');
    rho_perm = NaN(n_perm,size(y,1));
    for cur_perm = 1:n_perm
        xhat_perm = xhat(randperm(numel(xhat),numel(xhat)));
        rho_perm(cur_perm,:) = getRho(xhat_perm,y);
    end    
    rho_perm_all(:,cur_d) = max(abs(rho_perm));
    
end

%save off
fn = [savedir,rec_name,sprintf('_area%dmotif%d.mat',cur_a,cur_motif)]; 
codeversion = 'BehavioralNetworks_v2';
save(fn,'lag_all','rho_perm_all','alpha','rho_all','area_label','cur_motif','cur_a','codeversion','area_label');   
toc
    
end %function 
 

function [rho,lag] = getRho(xhat,y)
    %allows a +/- 1 second lag since we are just looking at broad
    %relationship to motoro activity and not a detailed temporal
    %relationship. 
    if size(y,1)>1
        [rho,lag] = arrayfun(@(n) (xcorr(xhat-nanmean(xhat),y(n,:)'-nanmean(y(n,:)),8,'normalized')), 1:size(y,1),'UniformOutput',0); 
        [rho,idx] = cellfun(@(x) max(x),rho);
        lag = arrayfun(@(x) lag{1}(x),idx);
    else   
        [rho,lag] = xcorr(xhat-nanmean(xhat),y-nanmean(y),7,'normalized');
        [rho,idx] = max(rho); 
        lag = lag(idx); 
    end
end
























