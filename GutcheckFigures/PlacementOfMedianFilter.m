%Place the median filter at each step in the processing pipeline to
%determine location that best resolves noise, but doesn't lose spatial
%information
%Generate plots at each point in the analysis pipeline
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SpeckleNoise';
if ~exist(savedir,'dir'); mkdir(savedir); end
%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_corticaldynamics';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string. 

fn_raw = GrabFiles('Pos0.ome.tif',0,{'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Mouse332_RestingState_NP_06_07_2021_1'});
fn_opts = GrabFiles('prepro_log',0,{'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Mouse332_RestingState_NP_06_07_2021_1'});
fn_raw = fn_raw{1};
opts = load(fn_opts{1});
opts = opts.prepro_log;


%% Preprocess
img_count = 251; 
stack = NaN(ceil(opts.crop_h/opts.spatial_bin_factor),...
    ceil(opts.crop_w/opts.spatial_bin_factor),...
    img_count);
stack_smooth = NaN(ceil(opts.crop_h/opts.spatial_bin_factor),...
    ceil(opts.crop_w/opts.spatial_bin_factor),...
    img_count);

%read tiff
warning ('off','all'); 
in_tiff = Tiff(fn_raw, 'r');
for i = 1:img_count %Crop and Allign Image
    in_tiff.setDirectory(i);   
    cur_img = in_tiff.read;   
    cur_img = imcrop(cur_img, opts.crop_position);
    if opts.angle <= 90 %if angled to right of yaxis
        cur_img = imrotate(cur_img,-opts.angle);
    else %if angled left of yaxis
        cur_img = imrotate(cur_img,180-opts.angle);
    end
    cur_img = padarray(cur_img,[60,60]); 
    cur_img = imcrop(cur_img,[opts.crop_cord(1),opts.crop_cord(2),...
        opts.crop_w,opts.crop_h]);

    %confirm correct dimensions since imcrop is inconsistent (see imcrop)
    cur_img = double(cur_img([1:opts.crop_h],[1:opts.crop_w])); 
    
    %Median filter
    cur_img_smooth = SpatialMedian(cur_img,5);
    
    %Mask and spatial bin
    cur_img_smooth = SpatialBin(cur_img_smooth,opts.spatial_bin_factor,opts.mask,1);
    cur_img = SpatialBin(cur_img,opts.spatial_bin_factor,opts.mask,1);

    %Remask to remove edge effects from mask
    cur_img_smooth(imresize(opts.mask,[68 68])==0)=0;
    cur_img(imresize(opts.mask,[68 68])==0)=0;
        
    stack(:,:,i) = cur_img;
    stack_smooth(:,:,i) = cur_img_smooth;
end
warning ('off','all'); 

data = {};

% Placed at very begining
%median filter in time
stack_smooth = SpatialMedian(stack,5);
stack_smooth = movmedian(stack_smooth,15,3);
dff = makeDFF(stack_smooth, opts); 
[fr,nanpxs] = deconvolveDFF(dff,gp);
% fr_binned = fr(1:2:end,:)+fr(2:2:end,:);
data{1} = conditionDffMat(fr,nanpxs);
temp = conditionDffMat(fr,nanpxs);
close all; for i = 1:96; imagesc(temp(:,:,i),[0 3]); title(sprintf('%d',i)); pause(0.05); end

% Placed after dff
dff = makeDFF(stack, opts); 
dff = SpatialMedian(dff);
[fr,nanpxs] = deconvolveDFF(dff,gp);
fr_binned = fr(1:2:end,:)+fr(2:2:end,:);
data{2} = conditionDffMat(fr_binned,nanpxs);

% Placed after deconvolution
dff = makeDFF(stack, opts); 
[fr,nanpxs] = deconvolveDFF(dff,gp);
fr = conditionDffMat(SpatialMedian(conditionDffMat(fr,nanpxs),gp.rawsmoothing_kernel(1)));
fr_binned = fr(1:2:end,:)+fr(2:2:end,:);
data{3} = conditionDffMat(fr_binned,nanpxs);


% Placed after binning (current version)
dff = makeDFF(stack, opts); 
[fr,nanpxs] = deconvolveDFF(dff,gp);
fr_binned = fr(1:2:end,:)+fr(2:2:end,:);
fr_binned = conditionDffMat(SpatialMedian(conditionDffMat(fr_binned,nanpxs),gp.rawsmoothing_kernel(1)));
data{4} = conditionDffMat(fr_binned,nanpxs);

% No filtering
dff = makeDFF(stack, opts); 
[fr,nanpxs] = deconvolveDFF(dff,gp);
fr_binned = fr(1:2:end,:)+fr(2:2:end,:);
data{5} = conditionDffMat(fr_binned,nanpxs);

%for comparison, compute the dff version
opts.method = 'movingavg';
dff = makeDFF(stack, opts); opts.method = 'zscore';
[fr,nanpxs] = deconvolveDFF(dff,gp);
fr_binned = fr(1:2:end,:)+fr(2:2:end,:);
dff_data = conditionDffMat(fr_binned,nanpxs);

%% For each version, plot example second, the residuals vs no filter, and the 5 components (spatial structure (PCA)) of the results and the residuals.

for i = 1:numel(data)
    cur_save_dir = [savedir filesep sprintf('version%d',i)];
    if ~exist(cur_save_dir,'dir'); mkdir(cur_save_dir); end
    close all
    %plot example second
    temp = data{i};
    cval = prctile(temp(:),99);
    for j = 1:15
        figure; hold on; imagesc(temp(:,:,j),[0 cval]); colormap magma; colorbar
        c = colorbar; ylabel(c, 'relative FR')
        set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
        axis equal;
        title(sprintf('Frame %d',j));
    end    
    saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('vers%d_signal',i),cur_save_dir,0); close all
    %plot the spatial components of the activity
    [coef, score, ~, ~, ~, mu] = pca(conditionDffMat(data{i})');
    score = conditionDffMat(score',nanpxs);
    for j = 1:5
        temp = score(:,:,j);
        figure; hold on; imagesc(temp,[prctile(temp(:),1),prctile(temp(:),99)]); colormap gray; colorbar
        set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
        axis equal;
        title(sprintf('Spatial component %d',j));
    end           
    saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('vers%d_spatialcomps',i),cur_save_dir,0); close all
    
    %plot the residuals
    if i <numel(data)
        resid = data{i}-data{end};
        cval = prctile(resid(:),99);
        for j = 1:15
            figure; hold on; imagesc(resid(:,:,j),[-cval cval]); colormap magma; colorbar
            c = colorbar; ylabel(c, 'residual relative FR')
            set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
            axis equal;
            title(sprintf('Residuals Frame %d',j));
        end
        saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('vers%d_resid',i),cur_save_dir,0); close all
        
        %plot the spatial components of the residuals
        [coef, score, ~, ~, ~, mu] = pca(conditionDffMat(resid)');
        score = conditionDffMat(score',nanpxs);
        for j = 1:5
            temp = score(:,:,j);
            figure; hold on; imagesc(temp,[prctile(temp(:),1),prctile(temp(:),99)]); colormap gray; colorbar
            set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
            axis equal;
            title(sprintf('Spatial component %d of residual',j));
        end    
        saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('vers%d_residcomps',i),cur_save_dir,0); close all
    end
          
end

%% plot the moving average for comparison. 
cur_save_dir = [savedir filesep 'medianintime'];
if ~exist(cur_save_dir,'dir'); mkdir(cur_save_dir); end
close all
%plot example second
temp = dff_data;
cval = prctile(temp(:),99);
for j = 1:15
    figure; hold on; imagesc(temp(:,:,j),[0 cval]); colormap magma; colorbar
    c = colorbar; ylabel(c, 'relative FR')
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('Frame %d',j));
end    
saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('signal',i),cur_save_dir,0); close all
%plot the spatial components of the activity
[coef, score, ~, ~, ~, mu] = pca(conditionDffMat(dff_data)');
score = conditionDffMat(score',nanpxs);
for j = 1:5
    temp = score(:,:,j);
    figure; hold on; imagesc(temp,[prctile(temp(:),1),prctile(temp(:),99)]); colormap gray; colorbar
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('Spatial component %d',j));
end           
saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('spatialcomps',i),cur_save_dir,0); close all



















