function [stack, opts] = CombineStacks(folder_path,parameter_class)
%use default gp.stack_suffix specified in general_params
gp = loadobj(feval(parameter_class));
[file_list,~] = GrabFiles(gp.stack_suffix,0,{folder_path});
num_files = numel(file_list); 

if num_files >1 %combine all the files in the folder
    %The first file has no identification number. So won't sort correctly. so remove and then add back
    indx = cellfun(@isempty, regexp(file_list,['(?<=_)\d+(?=',gp.stack_suffix,')'],'match','once'),'UniformOutput',0);
    first_file = file_list([indx{:}]==1);       
    file_list([indx{:}]==1)=[];%remove the first file from file_list (added back in a few lines)
    
    %catch errors
    if numel(first_file)~=1
        error('Too many or No files with no ID number. Check stack names');
    end

    %sort the remaining files according to their recording order
    [~, reindex] = sort(str2double(regexp(file_list, ['(?<=_)\d+(?=',gp.stack_suffix,')'], 'match', 'once' ))); %Sort by chunk 
    file_list = [first_file,file_list(reindex)];     %sort and add back first file

    %load all the files and compile
    stack = []; %no preallocation. whoops. 
    for cur_file = 1:num_files
       [~,name] = fileparts(file_list{cur_file});       
       fprintf('\n working on %s',name);
       if cur_file == 1 %load the first file
           temp = load(first_file{1});
           stack = cat(3,stack,temp.stack);
       else
           temp = load(file_list{cur_file});
           stack = cat(3,stack,temp.stack);
       end 
       if gp.delete_singlefiles %optionally delete the files
          delete(file_list{cur_file});
       end
    end %cur_file loop        

    opts = temp.opts; %go ahead and carry along the preprocesing opts 
else %just load and return the single file
    temp = load(file_list{1});
    stack = temp.stack; 
    opts = temp.opts; 
end


end

    
           
    
