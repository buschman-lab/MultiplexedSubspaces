function ImpactOfPCADenoise(stack)
%Camden MacDowell - timeless. Stack is a x,y,time imaging dff (e.g. post
%hemo correction). This looks at the impact on pca denoising

assert(size(stack,3)<=10000,'Use less timepoints. >10000 will take forever'); %keep stacks <10000 samples in length

gp = general_params;

data = stack; 
[x,y,z] = size(data);

%condition data and remove nan pxls
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

%level of denoising. In fraction pixels for gp.denoise_pxlfrac
lvl = [0, 0.01, 0.025, 0.05];
all_data = cell(1,numel(lvl));
for i = 1:numel(lvl)    
    fprintf('\n\t working on iter %d',i);

    %PCA Denoising
    data_temp = conditionDffMat(data,nanpxs,[], [x,y,z]);
    data_temp = DenoisePCA(data_temp,'denoise_pxlfrac',lvl(i));
    [data_temp,~] = conditionDffMat(data_temp);
    
    all_data{i} = data_temp;
end

%make figure;
figure('position',[223,558,1546,420]); hold on;
for i = 1:numel(lvl)    
    subplot(1,numel(lvl),i);
    if i==1
        imagesc(all_data{1}',[-0.25,0.25]); axis off
    else
        imagesc(all_data{i}'-all_data{1}',[-0.25,0.25]); axis off
    end
    title(sprintf('lvl %.2d',100*lvl(i)),'fontweight','normal');
end
sgtitle('comparison of decolved data with and without pca noise correction');


%make figure;
figure('position',[223,558,1546,420]); hold on;
for i = 1:numel(lvl)    
    subplot(1,numel(lvl),i);
    imagesc(all_data{i}',[0,0.5]); axis off
    title(sprintf('lvl %.2d',100*lvl(i)),'fontweight','normal');
end
sgtitle('comparison of decolved data with and without pca noise correction');

%compare explained variances
% pev = cellfun(@(X) CalculateExplainedVariance(data,data-X),all_data,'UniformOutput',1);
end
