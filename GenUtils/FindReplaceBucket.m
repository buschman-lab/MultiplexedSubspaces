function FindReplaceBucket(target_dir, filetypes)
%recursive search of files in target_dir to replace 'bucket' with 'cup'
%%INPUTS
%target_dir: full path of highest directory to search
%filetypes: 1xn cell array file type extensions ('m','py','txt')
%grabs all files of filetypes in target directory and all 
%subdirectories. Searches for 'bucket' and replaces with 'cup'
%recursive search for files
%tested on matlab 2020b on windows

file_list = cell(1,numel(filetypes)); 
file_list = cellfun(@(x) dir(fullfile(target_dir,['**\*.',x])),filetypes,'UniformOutput',0);%get list of files and folders in any subfolder
file_list = cat(1,file_list{:});
file_list = file_list(~[file_list.isdir]);  %remove any folders that made it into list (unlikely)

%loop through files and replace 
for cur_fn = 1:numel(file_list)
    % Read in the file as binary and convert to chars.
    fn = fullfile(file_list(cur_fn).folder,file_list(cur_fn).name);
    fid = fopen(fn);
    text = fread(fid, inf, '*char')';
    fclose(fid);
    % Find and replace.
    text = regexprep(text, 'bucket\', 'cup\');
    % Write out the new file.
    fid = fopen(fn, 'w');
    fwrite(fid, text);
    fclose(fid);
end

end