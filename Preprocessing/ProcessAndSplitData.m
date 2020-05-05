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
[~,~,z] = size(data);   
[data,nanpxs] = conditionDffMat(data); %nanpxs are the same on each iteration so fine to overwrite

%filter data
switch gp.w_deconvolution
    case 'simple_filter'
        fprintf('\n\tFiltering data')
        data = filterstack(data, opts.fps, gp.w_filter_freq, gp.w_filter_type, 1, 1);
        %Remove negative component of signal as in MacDowell 2020. Legacy. %Deconvolution with non-negative output is preferred (maintains more data, comes with own caveats). 
        data(data<0)=0;
    case 'lucric'
        fprintf('\n\tPerforming a Lucy-Goosey Deconvolution (Lucy-Richardson)\n')
        data = lucric(data,gp.d_gamma,gp.d_smooth,gp.d_kernel);
end

%normalize to 0 to 1 and make nonnegative
switch gp.w_normalization_method
    case 'pixelwise' %each between x and xth pixel intensity
        data_norm = NaN(size(data));
        for px = 1:size(data,2)
            data_norm(:,px) = (data(:,px)-prctile(data(:,px),gp.w_norm_val(1)))/(prctile(data(:,px),gp.w_norm_val(2))-prctile(data(:,px),gp.w_norm_val(1)));
        end             
    case 'full'
        data_norm = (data-prctile(data(:),gp.w_norm_val(1)))/(prctile(data(:),gp.w_norm_val(2))-prctile(data(:),gp.w_norm_val(1)));  
    case 'bounded'
        data_norm = (data)/(gp.w_norm_val(2)); %normalize between zero and the upper bound     
    otherwise
        error('Unknown normalization methods. Check general params')
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

end















