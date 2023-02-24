function mouse_num = MouseNumFromPath(folder_path,str_format)
%Camden MacDowell - timesless
%grabs mouse number from the folder path. Just add more switches under
%str_format to add more allowable types
%folder_path is cell array of paths
assert(iscell(folder_path),'Requires folder path to be cell input');

N=numel(folder_path);
mouse_num = ones(1,N);
for cur_fn = 1:N
    switch str_format
        case 'Mouse_' %search for Mouse and _ and get number in between
            [~,temp_str] = fileparts(folder_path{cur_fn});
            mouse_num(cur_fn) = str2double(erase(regexp(temp_str,'Mouse\d*_','match'),{'Mouse','_'}));
        case 'Mouse-' %search for Mouse and - and get number in between
            [~,temp_str] = fileparts(folder_path{cur_fn});
            mouse_num(cur_fn) = str2double(erase(regexp(temp_str,'Mouse\d*-','match'),{'Mouse','-'}));      
        otherwise
            error('unsupported str_format');
    end
end