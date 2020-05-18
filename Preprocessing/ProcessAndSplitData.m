function ProcessAndSplitData(fn,save_fn)
%Camden MacDowell - timeless
%filters, normalizes, and splits data in fn into training and test test
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
end
    
gp = general_params;

fprintf('\n\tLoading data')
%load data
temp = load(fn);
data = temp.stack;
opts = temp.opts;
clear temp;

%condition data and remove nan pxls
[x,y,z] = size(data);   
[data,nanpxs] = conditionDffMat(data); %nanpxs are the same on each iteration so fine to overwrite

%filter data
switch gp.w_deconvolution
    case 'simple_filter'
        fprintf('\n\tFiltering data')
        data = filterstack(data, opts.fps, gp.w_filter_freq, gp.w_filter_type, 1, 1);
        %Remove negative component of signal as in MacDowell 2020. Legacy. %Deconvolution with non-negative output is preferred (maintains more data, comes with own assumptions). 
        data(data<0)=0;
    case 'lucric'
        fprintf('\n\tPerforming a Lucy-Goosey Deconvolution (Lucy-Richardson)\n')
        data = lucric(data,gp.d_gamma,gp.d_smooth,gp.d_kernel);
end

%Denoise with PCA (removed banded pixels)
data = conditionDffMat(data,nanpxs,[], [x,y,z]);
data = DenoisePCA(data);
[data,~] = conditionDffMat(data);

%normalize to 0 to 1 
fprintf('\n\tPerforming %s normalization to %d value', gp.w_normalization_method, gp.w_norm_val);
switch gp.w_normalization_method
    case 'pixelwise' %each between x and xth pixel intensity
        data_norm = NaN(size(data));
        for px = 1:size(data,2)
            data_norm(:,px) = (data(:,px))/(prctile(data(:,px),gp.w_norm_val));
        end             
    case 'full' %normalize using the percentile of the maximum value per pixe. Use max for each pixel instead of all pixels because noise/artifact pixels may be continuously high throughout and warp distribution
        data_norm = data/prctile(max(data,[],1),gp.w_norm_val);          
    case 'bounded'
        data_norm = (data)/(gp.w_norm_val(2)); %normalize between zero and the upper bound     
    case 'none'
        data_norm = data;
    otherwise
        error('Unknown normalization method. Check general params')
end

%transpose (seqNMF operates rowwise)
data_norm = data_norm';

%Chunk (since MU is suboptimal for >4500 timepoints. As many as possible
num_chunks = floor(z/(gp.w_chunk_dur*opts.fps));
if ~isEven(num_chunks)% need even number
    num_chunks = num_chunks-1; 
end

%remove the remainder and reshape into chunks
data_trim = data_norm(:,1:end-mod(z,num_chunks*(gp.w_chunk_dur*opts.fps)));
data_trim = reshape(data_trim,[size(data_trim,1),num_chunks,(gp.w_chunk_dur*opts.fps)]);

%alternate testing and training
data_train = data_trim(:,1:2:num_chunks,:);
data_test = data_trim(:,2:2:num_chunks,:);

%save off the data in the scratch directory and the nanpxs
fprintf('\n\tSaving data')
save(save_fn,'data_norm','data_test','data_train','nanpxs','opts','gp','num_chunks','-v7.3')
fprintf('\n\tDONE')


%% for saveing off figure use; 
% [path, name] = fileparts(save_fn);
% saveCurFigs(gcf,'-dpng',sprintf('registration_%s',fn),[path filesep name, '_ProcessAndSplitFigures'],0); %close all;


end















