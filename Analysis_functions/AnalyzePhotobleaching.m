function AnalyzePhotobleaching(in_fns)

if nargin<1; [file_names, folder_names] = GrabFiles('_combined.mat'); end

data_mean = cell(1,numel(file_names));
data_variance = cell(1,numel(file_names));
for cur_file = 1:numel(file_names)
   temp = load(file_names{cur_file});
   temp = squeeze(nanmean(temp.stack,1:2))'; %get the mean trace
   data_mean{cur_file} = WindowedDist(temp,780,@(x) nanvar(x));
   data_variance{cur_file} = temp;
   
end

%pad to same size for each; 



end %function end

