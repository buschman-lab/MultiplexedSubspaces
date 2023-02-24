function data_all = LoadCorticalNetworks(reverseFlag)
%Camden 

if ispc,
    folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CorticalNetworks';
else
    folder = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CorticalNetworks';
end

%load by recording, adjust for missing areas, shift motif numbers
data_all = cell(1,6);
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    if reverseFlag
        [fn,~] = GrabFiles([rec_name,'REVERSE\w*.mat'],0,{folder});
    else
        [fn,~] = GrabFiles([rec_name,'area\w*.mat'],0,{folder});
    end
    data = cellfun(@(x) load(x),fn,'UniformOutput',1);
    
    %adjust motif numbers
    for i = 1:size(data,2)
        if data(i).cur_motif>1
            data(i).cur_motif = data(i).cur_motif-1;
        end
    end
    
    %adjust for missing areas in recs 3 and 4
    if cur_rec == 3 %RSP (4) and BFD (6)
        idx = [1,2,3,5,7,8];
        for i = 1:size(data,2)            
            data(i).cur_a = idx(data(i).cur_a);
        end
    elseif cur_rec ==4 %BFD (6)
        idx = [1,2,3,4,5,7,8];
        for i = 1:size(data,2)            
            data(i).cur_a = idx(data(i).cur_a);
        end
    end

    data_all{cur_rec} = data;
end