function data = LoadSubspaceData(flag)
fp = fig_params_cortdynamics;


data = cell(1,6);
switch flag
    case 'in_grouped'
        if ispc
            folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace';
        else
            folder = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CommunicationSubspace';
        end
        for cur_rec = 1:6
            rec_name = LoadDataDirectories(cur_rec);
            [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDmotif\d*.mat'],0,{folder}); 
            data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','inactive_idx','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif','st_depth'),fn);
        end
    case 'in_grouped_lag'
        folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace_Lagged';
        for cur_rec = 1:6
            rec_name = LoadDataDirectories(cur_rec);
            [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDmotif\d*.mat'],0,{folder}); 
            data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','inactive_idx','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif','st_depth'),fn);
        end   
    case 'out_grouped_lag'
        folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace_Lagged';
        for cur_rec = 1:6
            rec_name = LoadDataDirectories(cur_rec);
            [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDREVERSEmotif\d*.mat'],0,{folder}); 
            data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','inactive_idx','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif','st_depth'),fn);
        end           
%     case 'paired'
%         folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace';
%         for cur_rec = 1:6
%             rec_name = LoadDataDirectories(cur_rec);
%             [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_m\w*.mat'],0,{folder}); 
%             data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','inactive_idx','area_label','area_val','paired_areas','rrr_V','rrr_B','motif','st_depth'),fn);
%         end        
%     case 'out_grouped'
%         folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace';
%         for cur_rec = 1:6
%             rec_name = LoadDataDirectories(cur_rec);
%             [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDREVERSE\w*.mat'],0,{folder}); 
%             data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif','st_depth'),fn);
%         end
%     case 'in_grouped_meansubtract'
%         folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace_meansubtract';
%         for cur_rec = 1:6
%             rec_name = LoadDataDirectories(cur_rec);
%             [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDm\w*.mat'],0,{folder}); 
%             data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','inactive_idx','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif','st_depth'),fn);
%         end
%     case 'out_grouped_meansubtract'
%         folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace_meansubtract';
%         for cur_rec = 1:6
%             rec_name = LoadDataDirectories(cur_rec);
%             [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDREVERSE\w*.mat'],0,{folder}); 
%             data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif','st_depth'),fn);
%         end
%     case 'in_grouped_raw'
%         folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace';
%         for cur_rec = 1:6
%             rec_name = LoadDataDirectories(cur_rec);
%             [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDm\w*.mat'],0,{folder}); 
%             data{cur_rec} = cellfun(@(x) load(x,'st_norm','neu_area','motif'),fn(1)); %only need to load the first motifs, since all the same
%         end    
    otherwise
        error('unknown type of data to load');
end

%remove noise motif and null motif
for i = 1:6 
   idx = ismember(arrayfun(@(n) data{i}(n).motif, 1:size(data{i},2)),fp.noisemotif);
   data{i}(idx)=[];     
end

end %function end