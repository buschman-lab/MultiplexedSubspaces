function SpockCompareSubspaceBetas(shufseed,cur_rec,type)

%Camden macdowell - timeless
%quick code snippet to find jobs that failed
% [fn,~] = GrabFiles(['type\w*.mat'],0,{pwd});
% [~,fn] = cellfun(@(x) fileparts(x), fn, 'UniformOutput',0);
% %create true list
% temp = {};
% count=1;
% for type = 2:3
% for shuf = 1%:1000
% for rec = [1,2,3,4,5,6]
% a = sprintf('type%d_shuf%d_rec%d',type,shuf,rec);
% if sum(ismember(fn,a))==0
%     temp{count} = a;
%     count=count+1;
% end
% end
% end
% end


fprintf('working on seed %d',shufseed);

if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
end

%Load the data_rrr for our specific motifs from each rec
folder = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CommunicationSubspace';

rec_name = LoadDataDirectories(cur_rec);

switch type 
    case 1
        [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_m\w*.mat'],0,{folder}); 
        data = cellfun(@(x) load(x,'area_label','area_val','rrr_V','rrr_B','paired_areas','motif'),fn);
        %remove noise motif and null motif (if present)
        idx = ismember(arrayfun(@(n) data(n).motif, 1:size(data,2)),[2,16]);
        data(idx)=[];     

        fprintf('working on recording %d',cur_rec);
        [B,Bmatchdim] = SharedBetas(data,'B',10,shufseed);
        [V,Vmatchdim] = SharedBetas(data,'V',10,shufseed);
        [rsq,matchdim] = SharedBetas_paired(data,'projection',10,shufseed);
    case 2
        [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDm\w*.mat'],0,{folder}); 
        data = cellfun(@(x) load(x,'area_label','area_val','rrr_V','rrr_B','grouping','motif'),fn);        
        %remove noise motif and null motif (if present)
        idx = ismember(arrayfun(@(n) data(n).motif, 1:size(data,2)),[2,16]);
        data(idx)=[];     

        fprintf('working on recording %d',cur_rec);
        [B,Bmatchdim] = SharedBetas(data,'B',10,shufseed);
        [V,Vmatchdim] = SharedBetas(data,'V',10,shufseed);
        [rsq,matchdim] = SharedBetas(data,'projection',10,shufseed);

    case 3
        [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDR\w*.mat'],0,{folder}); 
        data = cellfun(@(x) load(x,'area_label','area_val','rrr_V','rrr_B','grouping','motif'),fn);        
        %remove noise motif and null motif (if present)
        idx = ismember(arrayfun(@(n) data(n).motif, 1:size(data,2)),[2,16]);
        data(idx)=[];     

        fprintf('working on recording %d',cur_rec);
        [B,Bmatchdim] = SharedBetas(data,'B',10,shufseed);
        [V,Vmatchdim] = SharedBetas(data,'V',10,shufseed);
        [rsq,matchdim] = SharedBetas(data,'projection_reverse',10,shufseed);

    otherwise
        error('unknown type');
end       


fn = sprintf('type%d_shuf%d_rec%d',type,shufseed,cur_rec);
% savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/SubspaceComparison_meansubtract/';
savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/SubspaceComparison/';
save([savedir,fn],'rsq','B','V','matchdim','Bmatchdim','Vmatchdim');


end






















