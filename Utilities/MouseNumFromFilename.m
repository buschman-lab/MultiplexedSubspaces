function mouse = MouseNumFromFilename(file_list)
    %Camden MacDowell - timeless
    %file names is a cell array of file paths where the name is
    %path/mousenum_otherinfo.mat;
    mouse = NaN(1,numel(file_list));
    for i = 1:numel(file_list)
        [~, temp] = fileparts(file_list{i});
        temp = regexp(temp,'_','split','once');
        mouse(i) = str2num(temp{1});
    end       
end
