function AnalyzePhotobleaching(in_fns)

if nargin<1; [file_names, folder_names] = GrabFiles('_combined.mat'); end

data = cell(1,numel(file_names));
for cur_file = 1:numel(file_names)
   temp = load(file_names{cur_file});
   temp = squeeze(nanmean(temp.stack,1:2))'; %get the mean trace
   data{cur_file} = WindowedDist(temp,780,@(x) nanvar(x));
end

%pad to same size for each; 



end %function end

