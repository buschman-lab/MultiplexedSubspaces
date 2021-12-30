%Generate plots at each point in the analysis pipeline
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ImpactOfPreprocessing';
%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_corticaldynamics';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string. 

fn_raw = GrabFiles('Pos0.ome.tif',0,{'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Mouse332_RestingState_NP_06_07_2021_1'});
fn_opts = GrabFiles('prepro_log',0,{'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Mouse332_RestingState_NP_06_07_2021_1'});
fn_raw = fn_raw{1};
opts = load(fn_opts{1});
opts = opts.prepro_log;

%plot the raw data 
warning ('off','all'); %Suppress missing tiff metadata warning
%get stack info (faster than tiff counter)
info = imfinfo(fn_raw); 
img_count = 30; %can't see much raw so just take a few
raw_stack = NaN(info(1).Height,info(1).Width,img_count);
%%
plot_window = 30:60; %plotting window. Can just choose a second or so at the begining. 

%read tiff
in_tiff = Tiff(fn_raw, 'r');
%Loop through each image
for cur_img_ind = plot_window
    in_tiff.setDirectory(cur_img_ind);   
    cur_img = in_tiff.read;   
    raw_stack(:,:,cur_img_ind) =  in_tiff.read;     
end %image loop
warning ('on','all');
% close all; for i = 1:img_count; imagesc(raw_stack(:,:,i)); title(sprintf('%d',i)); pause(); end 
close all; 
for i = 1:img_count
    figure; hold on; imagesc(raw_stack(:,:,i),[0 7000]); colormap gray; colorbar
    c = colorbar; ylabel(c, 'pxl intensity')
    set(gca,'ydir','reverse'); axis off; ylim([0 info(1).Height]); xlim([0 info(1).Width])
    axis equal;
    title(sprintf('Raw frame %d',i));
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('RawStack%d',i),savedir,0); close all

%plot the aligned and registered images
stack = PreProcess(fn_raw,opts);
close all; for i = 15:1000; imagesc(stack(:,:,i)); title(sprintf('%d',i)); pause(); end %offset to account for later deconvolution
close all; 
for i = plot_window
    figure; hold on; imagesc(stack(:,:,i),[0 2*10^5]); colormap gray; 
    c = colorbar; ylabel(c, 'pxl intensity')
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('Binned and Alligned frame %d',i));
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('Stack%d',i),savedir,0); close all

%plot the zscored data
dff = makeDFF(stack, opts); 
cval = prctile(dff(:),97);
close all; for i = 15:1000; imagesc(dff(:,:,i),[-1*cval cval]); title(sprintf('%d',i)); pause(0.05); end
close all; 
for i = plot_window
    figure; hold on; imagesc(dff(:,:,i),[-1*cval cval]); colormap gray; colorbar
    c = colorbar; ylabel(c, 'zscore')
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('DFS frame %d',i));
end
saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('dff%d',i),savedir,0); close all

%compare to traditional dff
opts.method = 'movingavg';
dff_mean = makeDFF(stack, opts);  opts.method = 'zscore';
cval = prctile(dff_mean(:),97);
close all; for i = 15:1000; imagesc(dff_mean(:,:,i),[-1*cval cval]); title(sprintf('%d',i)); pause(0.05); end
close all; 
for i = plot_window
    figure; hold on; imagesc(dff_mean(:,:,i),[-1*cval cval]); colormap gray; colorbar
    c = colorbar; ylabel(c, 'zscore')
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('DFS frame %d',i));
end
saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('dff_mean%d',i),savedir,0); close all

%plot the deconvolved data
[fr,nanpxs,x,y,~] = deconvolveDFF(dff,gp);
temp = conditionDffMat(fr,nanpxs);
close all; for i = 1:1000; imagesc(temp(:,:,i),[0 3.5]); title(sprintf('%d',i)); pause(0.05); end
for i = 1:30
    figure; hold on; imagesc(temp(:,:,i),[0 4]); colormap gray; colorbar
    c = colorbar; ylabel(c, 'relative FR')
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('FR frame %d',i));
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('fr%d',i),savedir,0); close all

%plot the binned data
fr_binned = fr(1:2:end,:)+fr(2:2:end,:);
opts.fps = floor(opts.fps/2);
z=size(fr_binned,1);
temp = conditionDffMat(fr_binned,nanpxs);
% close all; for i = 1:1000; imagesc(temp(:,:,i),[0 8]); title(sprintf('%d',i)); pause(0.1); end
for i = 1:15
    figure; hold on; imagesc(temp(:,:,i),[0 7]); colormap gray; colorbar
    c = colorbar; ylabel(c, 'relative FR')
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('FR_binned frame %d',i));
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('fr_binned%d',i),savedir,0); close all

%plot the median filter spatial smoothed data
fr_smooth = conditionDffMat(SpatialMedian(conditionDffMat(fr_binned,nanpxs),gp.rawsmoothing_kernel(1)));
temp = conditionDffMat(fr_smooth,nanpxs);
close all; for i = 1:1000; imagesc(temp(:,:,i),[0 8]); title(sprintf('%d',i)); pause(0.05); end
for i = 1:15
    figure; hold on; imagesc(temp(:,:,i),[0 7]); colormap gray; colorbar
    c = colorbar; ylabel(c, 'relative FR')
    set(gca,'ydir','reverse'); axis off; ylim([0 68]); xlim([0 68])
    axis equal;
    title(sprintf('FR_smoothed frame %d',i));
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('fr_smoothed%d',i),savedir,0); close all

%residuals
fr_resid = fr_smooth-fr_binned;
temp = conditionDffMat(fr_resid,nanpxs);
close all; for i = 1:1000; imagesc(temp(:,:,i),[-4 4]); title(sprintf('%d',i)); pause(0.05); end

[coef, score, ~, ~, ~, mu] = pca(fr_resid');
score = conditionDffMat(score',nanpxs);
close all; for i = 1:5; imagesc(score(:,:,i)); title(sprintf('%d',i)); pause(); end

[coef, score, ~, ~, ~, mu] = pca(fr_smooth');
score = conditionDffMat(score',nanpxs);
close all; for i = 1:5; imagesc(score(:,:,i)); title(sprintf('%d',i)); pause(); end



%%

