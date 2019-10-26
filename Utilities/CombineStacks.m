function [stack, opts] = CombineStacks(folder_path)
%use default gp.dff_suffix specified in general_params
gp = general_params; 
[file_list,~] = GrabFiles(gp.dff_suffix,0,{folder_path});
num_files = numel(file_list); 

if num_files >1 %combine all the files in the folder
    %The first file has no identification number. So won't sort correctly. 
    indx = cellfun(@isempty, regexp(file_list,['(?<=_)\d+(?=',gp.dff_suffix,')'], 'match', 'once'),'UniformOutput',0);
    first_file = file_list([indx{:}]==1);

    %catch errors
    if numel(first_file)~=1
        error('Too many or No files with no ID number. Check stack names');
    end

    %remove from file_list
    file_list([indx{:}]==1)=[];

    %sort the remaining files according to their recording order
    [~, reindex] = sort( str2double( regexp( file_list, ['(?<=_)\d+(?=',gp.dff_suffix,')'], 'match', 'once' ))); %Sort by block 
    file_list = file_list(reindex);

    %load all the files and compile
    stack = []; %no preallocation. whoops. 
    for cur_file = 1:num_files
       if cur_file == 1 %load the first file
           temp = load(first_file{1});
           stack = cat(3,stack,temp.stack);
       else
           temp = load(file_list{1});
           stack = cat(3,stack,temp.stack);
       end 
    end %cur_file loop 

    opts = temp.opts; %go ahead and carry along the preprocesing opts 
else %just load and return the single file
    temp = load(file_list{1});
    stack = temp.stack; 
    opts = temp.opts; 
end


end

    
           
    
