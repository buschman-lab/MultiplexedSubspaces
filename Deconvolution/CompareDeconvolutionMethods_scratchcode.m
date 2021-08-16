%plot correlation 
%todo: test the across animals
%test normalization methods
%switch figures to all median calculations
%plot the lag
%plot the shapes of the xcorr
%remove lr glm and narx and replace with an alternative 2 algs
%also compare max xcorrs


%% look at drift
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\full_deconvolutiondata.mat')
norm_method = 'std';
params.bindata = 0; %temporally bin the data?
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 600]; %depth from surface of probe
[~,st] = CompileData_deconvolution([],spike_opts_list,params);

%loop through each and get the magnitude of drift
b=[];
a=[];
for i = 1:numel(st)
    for j = 1:4
        temp = st{i}(30*60*60:end,j);        
        coef = polyfit(1:numel(temp),temp,1); 
        b(i,j) = coef(2);
        temp = xcorr(temp-nanmean(temp),600,'normalized');
        a(i,j) = temp(1);
    end
end




%%



%look at stats and switch to median 
%write the sophisticated figure code (treat each depth as a single figure)

%inspect raw calcium and decide if more preprocessing needed
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\full_deconvolutiondata.mat')
clear dff spike_opts_list

data = load(dff_list{1});
%pca to denoise
dff = data.dff;
stack_denoise = DenoisePCA(conditionDffMat(dff,nanpxs));
[coef,score,~,~,explained,mu]=pca(dff);
%take first 500PCs
dff_denoised = score(:,1:1000)*coef(:,1:1000)' + repmat(mu,size(score,1),1);

nanpxs = data.nanpxs;

%% manually inspect the data
figure;
for i = 2:6
%     temp = conditionDffMat(data{i}.dff,data{i}.nanpxs);
%     temp = conditionDffMat(dff_denoised,nanpxs);
%     temp = conditionDffMat(dff,nanpxs);
    temp = stack_denoise;
    temp = temp(:,:,1:2:end)+temp(:,:,2:2:end);
    for j = 80000:81000 %grab from the start, middle, and end of the recording
       imagesc(temp(:,:,j),[-5 5]); colormap magma; 
       title(sprintf('rec %d frame %d',i,j));
       pause(0.01)
    end
end 
%plot range over recording time, plot the correlational structure of
%different seed areas


%% scratch code for plotting the deconvolve comparisons with nonparameteric !!!!!!!!!!!!!!!!!!!!!!!!!!
% fn = GrabFiles('within_compare\w*mapminmax\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn = GrabFiles('within_compare\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
% fn = GrabFiles('within_compare\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\timebinned'});
rng('default') % For reproducibility
%loop through depths and get the average correlation
rho_avg = cell(1,numel(fn));
rho_err = cell(1,numel(fn));
mse_avg = cell(1,numel(fn));
mse_err = cell(1,numel(fn));
sr = nan(2,numel(fn));
t = nan(2,numel(fn));
for i = 1:numel(fn)
    data = load(fn{i});
    deconv_type = {'narx','lr_gcamp','lr_glm','glm','feedforward','none'}; 
    deconv_type = repmat(deconv_type,size(data.test_idx,1),1);
    %remove skipped 
    data.test_idx(data.bad_probe_idx,:)=[];
    data.train_idx(data.bad_probe_idx,:)=[];
    deconv_type(data.bad_probe_idx,:)=[];
    deconv_type = cat(1,deconv_type(:));
    %remove all occurances
    type = unique(deconv_type,'stable');   
    rho_avg{i} = NaN(1,numel(type));
    rho_err{i} = NaN(1,numel(type));
    mse_avg{i} = NaN(1,numel(type));
    mse_err{i} = NaN(1,numel(type));    
    for j = 1:numel(type)
        idx = strcmp(deconv_type,type{j});
        rho_avg{i}(j) = nanmean([data.deconv_stats(idx).rho]);
        rho_err{i}(j) = sem([data.deconv_stats(idx).rho]');
        mse_avg{i}(j) = nanmean([data.deconv_stats(idx).err]);
        mse_err{i}(j) = sem([data.deconv_stats(idx).err]');
    end       
    sr(1,i) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
    sr(2,i) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
    [~,t(1,i)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
    [~,t(2,i)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
end

rho_avg = cat(1,rho_avg{:});
rho_err = cat(1,rho_err{:});
mse_avg = cat(1,mse_avg{:});
mse_err = cat(1,mse_err{:});

%across all depths plot the correlation by depth
figure; hold on; 
for i = 1:size(rho_avg,1)
    depth = data.depths;
    x = 1:size(rho_avg,2);    
    b=bar(x+(8*(i-1)),rho_avg(i,:),'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
    b.CData = getColorPalet(6);
    er = errorbar(x+(8*(i-1)),rho_avg(i,:),rho_err(i,:),rho_err(i,:),'linestyle','none','linewidth',2,'color','k');    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';      
end
set(gca,'xtick',[3.5:8:i*8],'xticklabels',depth(:,2))
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('Rho')

%across all depths plot the correlation by depth
figure; hold on; 
for i = 1:size(mse_avg,1)
    depth = data.depths;
    x = 1:size(mse_avg,2);    
    b=bar(x+(8*(i-1)),mse_avg(i,:),'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
    b.CData = getColorPalet(6);
    er = errorbar(x+(8*(i-1)),mse_avg(i,:),mse_err(i,:),mse_err(i,:),'linestyle','none','linewidth',2,'color','k');     
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';      
end
set(gca,'xtick',[3.5:8:i*8],'xticklabels',depth(:,2))
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('MSE')



%% Compare within animal, across sites
% data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\within_xsitesmapminmax.mat');
% data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\within_xsitesstd.mat');
data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\timebinned\within_xsitesstd.mat');
deconv_type = {'narx','lr_gcamp','lr_glm','glm','feedforward','none'}; 
deconv_type = repmat(deconv_type,size(data.test_idx,1),1);
%remove skipped 
data.test_idx(data.bad_probe_idx,:)=[];
data.train_idx(data.bad_probe_idx,:)=[];
deconv_type(data.bad_probe_idx,:)=[];
deconv_type = cat(1,deconv_type(:));
%remove all occurances
type = unique(deconv_type,'stable');   
rho_avg = NaN(1,numel(type));
rho_err = NaN(1,numel(type));
mse_avg = NaN(1,numel(type));
mse_err = NaN(1,numel(type));
lag_avg = NaN(1,numel(type));
lag_sem = NaN(1,numel(type));
skew_avg = NaN(1,numel(type));
skew_sem = NaN(1,numel(type));
xplot = cell(1,numel(type));
for j = 1:numel(type)
    idx = strcmp(deconv_type,type{j});
    rho_avg(j) = nanmean([data.deconv_stats(idx).rho]);
    rho_err(j) = sem([data.deconv_stats(idx).rho]');
    mse_avg(j) = nanmean([data.deconv_stats(idx).err]);
    mse_err(j) = sem([data.deconv_stats(idx).err]');    
    lag_avg(j) = nanmean([data.deconv_stats(idx).x_lag]);
    lag_sem(j) = sem([data.deconv_stats(idx).x_lag]');    
    skew_avg(j) = nanmean([data.deconv_stats(idx).s]);
    skew_sem(j) = sem([data.deconv_stats(idx).s]'); 
    xplot{j} = [data.deconv_stats(idx).xcorrvect]';
end       
clear sr t
sr(1) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
sr(2) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
[~,t(1)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
[~,t(2)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
sr(3) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).s]',[data.deconv_stats(strcmp(deconv_type,type{5})).s]');
[~,t(3)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).s]',[data.deconv_stats(strcmp(deconv_type,type{5})).s]');

%per type, get the plot of the cross corr
figure; hold on; 
col =  getColorPalet(6);
p = cell(1,6);
for i = 1:numel(xplot)
    x = -60:60;
    y = nanmean(xplot{i});
    yerr = sem(xplot{i});
    shadedErrorBar(x,y,yerr,'lineprops',{'linestyle','-','color',col(i,:)});
    p{i} = plot(x,y,'linestyle','-','color',col(i,:));
end
legend([p{:}],type)
title('xcorr across sites')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(rho_avg,2);    
b=bar(x,rho_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,rho_avg,rho_err,rho_err,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('Rho')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(mse_avg,2);    
b=bar(x,mse_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,mse_avg,mse_err,mse_err,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('MSE')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(lag_avg,2);    
b=bar(x,lag_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,lag_avg,lag_sem,lag_sem,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('lag')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(skew_avg,2);    
b=bar(x,skew_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,skew_avg,skew_sem,skew_sem,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('delta Skew')

%need to add comparisons across specific brain areas
%also in stats and the similarity in the distribution shape and lag 

%% Compare across recordings in same animal
% data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\xrec_xsitesmapminmax.mat');
% data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\xrec_xsitesstd.mat');
data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\timebinned\xrec_xsitesstd.mat');
deconv_type = {'narx','lr_gcamp','lr_glm','glm','feedforward','none'}; 
deconv_type = repmat(deconv_type,size(data.test_idx,1),1);

%remove skipped 
data.test_idx(data.bad_probe_idx,:)=[];
data.train_idx(data.bad_probe_idx,:)=[];
deconv_type(data.bad_probe_idx,:)=[];
deconv_type = cat(1,deconv_type(:));
%remove all occurances
type = unique(deconv_type,'stable');   
rho_avg = NaN(1,numel(type));
rho_err = NaN(1,numel(type));
mse_avg = NaN(1,numel(type));
mse_err = NaN(1,numel(type));
lag_avg = NaN(1,numel(type));
lag_sem = NaN(1,numel(type));
skew_avg = NaN(1,numel(type));
skew_sem = NaN(1,numel(type));
xplot = cell(1,numel(type));
for j = 1:numel(type)
    idx = strcmp(deconv_type,type{j});
    rho_avg(j) = nanmean([data.deconv_stats(idx).rho]);
    rho_err(j) = sem([data.deconv_stats(idx).rho]');
    mse_avg(j) = nanmean([data.deconv_stats(idx).err]);
    mse_err(j) = sem([data.deconv_stats(idx).err]');    
    lag_avg(j) = nanmean([data.deconv_stats(idx).x_lag]);
    lag_sem(j) = sem([data.deconv_stats(idx).x_lag]');    
    skew_avg(j) = nanmean([data.deconv_stats(idx).s]);
    skew_sem(j) = sem([data.deconv_stats(idx).s]'); 
    xplot{j} = [data.deconv_stats(idx).xcorrvect]';
end       
clear sr t
sr(1) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
sr(2) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
[~,t(1)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
[~,t(2)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
sr(3) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).s]',[data.deconv_stats(strcmp(deconv_type,type{5})).s]');
[~,t(3)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).s]',[data.deconv_stats(strcmp(deconv_type,type{5})).s]');

%per type, get the plot of the cross corr
figure; hold on; 
col =  getColorPalet(6);
p = cell(1,6);
for i = 1:numel(xplot)
    x = -60:60;
    y = nanmean(xplot{i});
    yerr = sem(xplot{i});
    shadedErrorBar(x,y,yerr,'lineprops',{'linestyle','-','color',col(i,:)});
    p{i} = plot(x,y,'linestyle','-','color',col(i,:));
end
legend([p{:}],type)
title('xcorr across recs')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(rho_avg,2);    
b=bar(x,rho_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,rho_avg,rho_err,rho_err,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('Rho')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(mse_avg,2);    
b=bar(x,mse_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,mse_avg,mse_err,mse_err,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('MSE')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(lag_avg,2);    
b=bar(x,lag_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,lag_avg,lag_sem,lag_sem,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('lag')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(skew_avg,2);    
b=bar(x,skew_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,skew_avg,skew_sem,skew_sem,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('delta Skew')

%% Compare across animals, same 'site' =
% data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\xanimals_withinsitesmapminmax.mat');
% data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\xanimals_withinsitesstd.mat');
data = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\timebinned\xanimals_withinsitesstd.mat');
deconv_type = {'narx','lr_gcamp','lr_glm','glm','feedforward','none'}; 
deconv_type = repmat(deconv_type,size(data.test_idx,1),1);

%remove skipped 
data.test_idx(data.bad_probe_idx,:)=[];
data.train_idx(data.bad_probe_idx,:)=[];
deconv_type(data.bad_probe_idx,:)=[];
deconv_type = cat(1,deconv_type(:));
%remove all occurances
type = unique(deconv_type,'stable');   
rho_avg = NaN(1,numel(type));
rho_err = NaN(1,numel(type));
mse_avg = NaN(1,numel(type));
mse_err = NaN(1,numel(type));
lag_avg = NaN(1,numel(type));
lag_sem = NaN(1,numel(type));
skew_avg = NaN(1,numel(type));
skew_sem = NaN(1,numel(type));
xplot = cell(1,numel(type));
for j = 1:numel(type)
    idx = strcmp(deconv_type,type{j});
    rho_avg(j) = nanmean([data.deconv_stats(idx).rho]);
    rho_err(j) = sem([data.deconv_stats(idx).rho]');
    mse_avg(j) = nanmean([data.deconv_stats(idx).err]);
    mse_err(j) = sem([data.deconv_stats(idx).err]');    
    lag_avg(j) = nanmean([data.deconv_stats(idx).x_lag]);
    lag_sem(j) = sem([data.deconv_stats(idx).x_lag]');    
    skew_avg(j) = nanmean([data.deconv_stats(idx).s]);
    skew_sem(j) = sem([data.deconv_stats(idx).s]'); 
    xplot{j} = [data.deconv_stats(idx).xcorrvect]';
end       
clear sr t
sr(1) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
sr(2) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
[~,t(1)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).rho]',[data.deconv_stats(strcmp(deconv_type,type{5})).rho]');
[~,t(2)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).err]',[data.deconv_stats(strcmp(deconv_type,type{5})).err]');
sr(3) = signrank([data.deconv_stats(strcmp(deconv_type,type{4})).s]',[data.deconv_stats(strcmp(deconv_type,type{5})).s]');
[~,t(3)] = ttest([data.deconv_stats(strcmp(deconv_type,type{4})).s]',[data.deconv_stats(strcmp(deconv_type,type{5})).s]');

%per type, get the plot of the cross corr
figure; hold on; 
col =  getColorPalet(6);
p = cell(1,6);
for i = 1:numel(xplot)
    x = -60:60;
    y = nanmean(xplot{i});
    yerr = sem(xplot{i});
    shadedErrorBar(x,y,yerr,'lineprops',{'linestyle','-','color',col(i,:)});
    p{i} = plot(x,y,'linestyle','-','color',col(i,:));
end
legend([p{:}],type)
title('xcorr across animals')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(rho_avg,2);    
b=bar(x,rho_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,rho_avg,rho_err,rho_err,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('Rho')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(mse_avg,2);    
b=bar(x,mse_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,mse_avg,mse_err,mse_err,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('MSE')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(lag_avg,2);    
b=bar(x,lag_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,lag_avg,lag_sem,lag_sem,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('lag')

%across all depths plot the correlation by depth
figure; hold on; 
x = 1:size(skew_avg,2);    
b=bar(x,skew_avg,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor',[0.25 0.25 0.25]);
b.CData = getColorPalet(6);
er = errorbar(x,skew_avg,skew_sem,skew_sem,'linestyle','none','linewidth',2,'color','k');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';      
set(gca,'xtick',[1:6],'xticklabels',type)
title(sprintf('%s ',type{:}),'fontweight','normal','Interpreter','none');
xlabel('Depth (um)')
ylabel('delta Skew')

%% To do. 
%1) Check out the glm deconvolution
%2) Figure out the best normaliation (maybe follow what they did?)
%2) Decide on best number of pixels to use (rerun spatial ephys correlation
%maps)
%3) Trouble shoot other comparisons
%4) Make prelim figures of the other comparisons
%5) Look at the raw imaging data and really see how it looks/do any
%filtering or anything needed
%6) Potential rerun and redo of the prelim figures
%7) Present to Tim and settle on figures
%8) Sketch out figures and begin to write


%%
temp = cat(1,trained_opts{:});
temp = cat(2,temp(:,4).glmkernel);
figure; hold on; 
for i = 1:size(temp,2)
   a = temp(:,i);
%    a = (a-min(a))/(max(a)-min(a));
   a = a/std(a);
   temp(:,i)=a;
   plot(a,'color',[0.75 0.75 0.75],'linewidth',1)   
end
plot(nanmean(temp,2),'color','k','linewidth',2);






%% troubleshooting things
%split train/test (additional train,validation,test occurs within train)


load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\full_deconvolutiondata.mat')
spike_opts_list = spike_opts_list(4);
dff_list = dff_list(4);

%compile imaging data
fprintf('\n\t Compiling Imaging Data')
params.bindata = 0; %temporally bin the data?
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %m

%run for a given block (depth). Takes about 4-12 hours per depth so do this way
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 600]; %depth from surface of probe
%compile spiking data
dff = dff(4);
[~,st] = CompileData_deconvolution([],spike_opts_list,params);
n = floor(size(dff{1},1)*3/4);
dff_train = cellfun(@(x) mapstd(x(1:n,:)')',dff,'UniformOutput',0);
dff_test = cellfun(@(x) mapstd(x(n+1:end,:)')',dff,'UniformOutput',0);            
st_train = cellfun(@(x) mapstd(x(1:n,:)')',st,'UniformOutput',0);
st_test = cellfun(@(x) mapstd(x(n+1:end,:)')',st,'UniformOutput',0);  
trained_opts = Deconvolve_Train(dff_train,st_train,'glm',100,params.bindata);

n_rec=6;
cur_rec = 1;
cur_probe = 1;
trained_opts_cur_rec = struct();
trace = dff_train{cur_rec};
%     spikes = st_train{cur_rec}/std(st_train{cur_rec});   
%     spikes_test = st_test{cur_rec}/std(st_test{cur_rec}); 
spikes = st_train{cur_rec};   
spikes_test = st_test{cur_rec}; 
trace_test = dff_test{cur_rec};
n_probe = size(trace,2);

trace_probe = trace(:,cur_probe); 
spikes_probe = spikes(:,cur_probe);
spikes_probe_test = spikes_test(:,cur_probe);
trace_probe_test = trace_test(:,cur_probe);
gamma = [0.60:0.005:1];
rho = NaN(1,numel(gamma));
e = NaN(1,numel(gamma));
for i = 1:numel(gamma)
    stPred = lucric(trace_probe-min(trace_probe),gamma(i),1,60); 
    stats = deconvolutionfitStats(stPred,spikes_probe);
    rho(i) = stats.rho;
end
figure; hold on; 
plot(rho,'linewidth',2,'color',[0.5 0.5 0.5]);
[a,b] = max(rho)
gamma(b)
ylabel('rho'); xlabel('gamma'); set(gca,'XTick',1:20:numel(gamma),'XTickLabel',gamma(1:20:end))
xlim([1 numel(gamma)])
hold on; plot(b,a,'marker','o','color','k','linewidth',2)
idx = find(gamma==0.96);
hold on; plot(idx,rho(idx),'marker','o','color','r','linewidth',2)
title('optimizing lucrid gamma')

%         trace_probe = mapminmax(trace(:,cur_probe)',0,1)'; 
%         spikes_probe = mapminmax(spikes(:,cur_probe)',0,1)';
%         spikes_probe_test = mapminmax(spikes_test(:,cur_probe)',0,1)';
%         trace_probe_test = mapminmax(trace_test(:,cur_probe)',0,1)';        
%         
%         trace_probe = trace(:,cur_probe)/std(trace(:,cur_probe)); 
%         spikes_probe = spikes(:,cur_probe)/std(spikes(:,cur_probe));
%         spikes_probe_test = spikes_test(:,cur_probe)/std(spikes_test(:,cur_probe));
%         trace_probe_test = trace_test(:,cur_probe)/std(trace_test(:,cur_probe));        
%glm fit
win = 61; %size of the kernel window in frames (def = 1 sec=30)
predictors = createRollingWindow(trace_probe, win); %t-n:t-1
response =  spikes_probe(ceil(win/2):end-floor(win/2)); % get the middle timepoint in window  
%         response = response/std(response);
kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
glmkernel = kernel(2:end); %remove intercept  
inter=kernel(1);

stPred = convn(padarray(trace_probe',[0,floor(length(glmkernel)/2)],'replicate','both')',glmkernel,'valid')+inter; 
stPred = convn(padarray(trace_probe',[0,floor(length(glmkernel)/2)],'replicate','both')',flipud(glmkernel),'valid');  
% stPred(1)=[];
%test on the training data
%if crappy, see how improves with different convolution... that could be
%the problem
stats = deconvolutionfitStats(stPred,spikes_probe);

stPred = convn(padarray(trace_probe_test',[0,floor(length(glmkernel)/2)],'replicate','both')',flipud(glmkernel),'valid');  
% stPred(1)=[];
%test on the training data
%if crappy, see how improves with different convolution... that could be
%the problem
stats = deconvolutionfitStats(stPred,spikes_probe_test);

%also check out the LR deconv for same reason

kernel = glmkernel; 
kernel = kernel-min(kernel); %requires positive
stPred = deconvlucy(trace_probe_test-min(trace_probe_test),kernel); 
stats = deconvolutionfitStats(stPred,spikes_probe_test);

rng('default');
params = [];
params.n_hiddenlayer = 10; %neurons in hidden layer
params.trainFcn = 'trainlm'; % 'trainlm' works best
params.win = 60; %size of timepoints to use. 60 is best
params.verbose = 0; 
params.hiddenfnc = 'tansig'; %tansig, softmax, and radbasn, are all good
params.outputfnc = 'purelin'; %pure linear performs best, poslin is comparable but forces postive outpu               
net = train_feedforward_nn(trace_probe',spikes_probe',params); %train 

x = createRollingWindow(trace_probe', params.win)'; %t-n:t-1        
stPred = net(x)';        
% timepoints(ceil(params.win/2):end-floor(params.win/2))=1; % get the middle timepoint in window 
stats = deconvolutionfitStats(stPred,spikes_probe(ceil(params.win/2):end-floor(params.win/2)))

x = createRollingWindow(trace_probe_test', params.win)'; %t-n:t-1        
stPred = net(x)';        
% timepoints(ceil(params.win/2):end-floor(params.win/2))=1; % get the middle timepoint in window 
stats = deconvolutionfitStats(stPred,spikes_probe_test(ceil(params.win/2):end-floor(params.win/2)))




%%


%load file names
%Todo. Another network to test would be train on t-30:t+30 (no
%autoregression)
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\temp_data\decondata.mat')
%grab recordings
dff_list = GrabFiles('dff_combined.mat',1,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging'});
spike_opts_list = GrabFiles('ap_opts.mat',1,{'H:\'});

%% Load example data
% load the data
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 250]; %depth from surface of probe
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
params.splittype = 'half'; %type of training/test split to use (current just first and second half of recording
[dff,st,~] = CompileData_deconvolution(dff_list(4),spike_opts_list(4),params);

%% Try out a feedforward network
rng('default');
% best parameters for the feedforward
params.n_hiddenlayer = 10; %neurons in hidden layer
params.trainFcn = 'trainlm'; % 'trainlm' works best
params.win = 60; %size of timepoints to use. 60 is best
params.verbose = 0; 
params.hiddenfnc = 'tansig'; %tansig, softmax, and radbasn, are all good
params.outputfnc = 'purelin'; %pure linear performs best, poslin is comparable but forces postive output
cur_probe = 3;
n = size(dff{1},1)*3;
trace_probe_train = dff{1}(1:floor(n/4),cur_probe); 
spikes_probe_train = st{1}(1:floor(n/4),cur_probe); 
trace_probe_test = dff{1}(floor(n/4)+1:end,cur_probe); 
spikes_probe_test = st{1}(floor(n/4)+1:end,cur_probe); 
[netc,stats] = train_feedforward_nn(trace_probe_train',spikes_probe_train',params); %train 
%test it
x = createRollingWindow(trace_probe_train', 60)'; %t-n:t-1
t =  spikes_probe_train(ceil(60/2):end-floor(60/2))'; % get the middle timepoint in window  

y_train = netc(x);
figure; hold on; 
plot(t); plot(y_train)
figure; plot(t,y_train,'.')
rho = corr(t',y_train')


x = createRollingWindow(trace_probe_test', 60)'; %t-n:t-1
t =  spikes_probe_test(ceil(60/2):end-floor(60/2))'; % get the middle timepoint in window  
y_train = netc(x);
figure; hold on; 
plot(t); plot(y_train)
figure; plot(t,y_train,'.')
rho = corr(t',y_train')


%% Choose best network
comp_fn = GrabFiles('block\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\temp_data'});
data = cellfun(@(x) load(x,'rho','e','params','block','rho_gen','e_gen'),comp_fn,'UniformOutput',0);

[block,idx] = sort(cellfun(@(x) x.block, data,'UniformOutput',1),'ascend');
rho = cellfun(@(x) x.rho, data(idx),'UniformOutput',1);
rho_gen = cellfun(@(x) x.rho_gen, data(idx),'UniformOutput',1);
e = cellfun(@(x) x.e, data(idx),'UniformOutput',1);
e_gen = cellfun(@(x) x.e_gen, data(idx),'UniformOutput',1);
params = cellfun(@(x) x.params, data(idx),'UniformOutput',0);
figure; plot(rho,'k'); hold on; plot(rho_gen,'r')
figure; hold on; plot(log10(e),'k'); plot(log10(e_gen),'r');
[a,b] = maxk(rho,10);
[a,b] = maxk(rho_gen,10);
[ea,eb] = mink(e_gen,10);
temp = cat(1,params{eb});

%best parameters = 



%%
X = tonndata(trace_probe_train',true,false);
rho = NaN(1,numel(func_list));
e = NaN(1,numel(func_list));
func_list = {'compet','elliotsig','hardlim','hardlims','logsig','netinv','poslin','purepin','radbas','radbasn','satlin','softmax','tansig','tribas'};
for cur_func = 1:numel(func_list)
    try 
        tic
        rng('default');
        fprintf('\n\t working on function %d of %d \n',cur_func, numel(func_list));
        params.outputfnc  = func_list{cur_func};
        [~,~,netc,~,~] = train_narx_nn(trace_probe_train',spikes_probe_train',params); %train   
        [x,xic,~,~] = preparets(netc,X,{});
        aic = cat(1,repmat({zeros(params.n_hiddenlayer,1)},1,params.feedbackdelay),repmat({zeros(1,1)},1,params.feedbackdelay)); %init feedback condition
        y_train = netc(x,xic,aic);        
        y_train = ([y_train{:}]');
        temp = spikes_probe_train(params.inputdelay+1:end);
        rho(cur_func) = corr(y_train,temp(1:numel(y_train)));
        e(cur_func) = mse(y_train,T);    
        fprintf('\n\t Function %s Rho%0.2g took %0.0f secs \n',func_list{cur_func},rho(cur_func), toc);
    catch
        fprintf('no luck with %s',func_list{cur_func});
    end        
end
%%
figure; hold on; 
plot(rho); 
e_train = e; 



%train network 
% [~,train_stats,netc,xic,aic] = train_narx_nn(flipud(trace_probe_train)',flipud(spikes_probe_train)',params); %train

% %     compet - Competitive transfer function.
% %     elliotsig - Elliot sigmoid transfer function.
% %     hardlim - Positive hard limit transfer function.
% %     hardlims - Symmetric hard limit transfer function.
% %     logsig - Logarithmic sigmoid transfer function.
% %     netinv - Inverse transfer function.
% %     poslin - Positive linear transfer function.
% %     purelin - Linear transfer function.
% %     radbas - Radial basis transfer function.
% %     radbasn - Radial basis normalized transfer function.
% %     satlin - Positive saturating linear transfer function.
% %     satlins - Symmetric saturating linear transfer function.
% %     softmax - Soft max transfer function.
% %     tansig - Symmetric sigmoid transfer function.
% %     tribas - Triangular basis transfer function.




%% Example figures using recording 4, probe 1, 1:250uM depth 
% train the network
rng('default');
params.n_hiddenlayer = 20; %neurons in hidden layer
params.trainFcn = 'trainlm'; % 'trainbr'
params.inputdelay = 30; %number of timepoints for prediction
params.feedbackdelay = 5; %number of timepoints for predictionparams.verbose = 0; 
params.hiddenfnc = 'tansig';
params.outputfnc = 'purelin';
params.verbose = 0; 
cur_probe = 1; 

%split data 
n = size(dff{1},1)*3;
trace_probe_train = dff{1}(1:floor(n/4),cur_probe); 
spikes_probe_train = st{1}(1:floor(n/4),cur_probe); 
trace_probe_test = dff{1}(floor(n/4)+1:end,cur_probe); 
spikes_probe_test = st{1}(floor(n/4)+1:end,cur_probe); 

%train network 
[~,train_stats,netc,xic,aic] = train_narx_nn((trace_probe_train)',(spikes_probe_train)',params); %train
% [~,train_stats,netc,xic,aic] = train_narx_nn(flipud(trace_probe_train)',flipud(spikes_probe_train)',params); %train

%closed loop network performance
%train data
% X = tonndata((trace_probe_train([train_stats.train_indx{:}]==1))',true,false);
% T = tonndata((spikes_probe_train([train_stats.train_indx{:}]==1))',true,false);
X = tonndata(trace_probe_train',true,false);
T = tonndata(spikes_probe_train',true,false);
[x,xic,~,t_train] = preparets(netc,X,{},T);
aic = cat(1,repmat({zeros(params.n_hiddenlayer,1)},1,params.feedbackdelay),repmat({zeros(1,1)},1,params.feedbackdelay)); %init feedback condition
y_train = netc(x,xic,aic);   
x_train = ([x{:}]');
y_train = ([y_train{:}]');
t_train = ([t_train{:}]');


%test data
X = tonndata(flipud(trace_probe_test)',true,false);
T = tonndata(flipud(spikes_probe_test)',true,false);
[x,~,~,t_test] = preparets(netc,X,{},T);
y_test = netc(x,xic,aic);    
x_test = flipud([x{:}]');
y_test = flipud([y_test{:}]');
t_test = flipud([t_test{:}]');

%compare with other methods on testing data
win = 60; %size of the kernel window in frames (def = 1 sec=30)
predictors = createRollingWindow(trace_probe_test,win); %t-n:t-1
response =  spikes_probe_test(ceil(win/2):end-floor(win/2)); % get the middle timepoint in window  
kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
mu = kernel(1);
kernel = kernel(2:end);
y_glm = convn(trace_probe_test,kernel,'same')+mu;
y_glm_lr = deconvlucy(trace_probe_test-min(trace_probe_test), kernel-min(kernel))+mu; 
y_lr = lucric(trace_probe_test-min(trace_probe_test),0.95,1,30)+mu;

%Figures
%plot middle 1 minute chunck
t = [10*30*20+1:20*30*11]; %min 10 to 11
x = [1:numel(t)]/30; %in seconds
figure('units','centimeters','position',[13 13 15 4]); 
hold on; plot(x,t_train(t),'color',[0.25 0.25 0.25,0.25],'linestyle','-','linewidth',1); 
hold on; plot(x,y_train(t),'color',[0.612 0.25 0.45],'linestyle','-','linewidth',1);
ylabel('Firing Rate'); 
yyaxis right
hold on; plot(x,x_train(t),'color',[0 0.25 0.612 0.25],'linestyle','-','linewidth',1);
ylabel('DFF','Interpreter','latex'); 
set(gca,'YColor','k');
xlabel('Time (s)');
title('Exampled Data from Training Set');

%testing data
figure('units','centimeters','position',[13 13 15 4]); 
hold on; plot(x,t_test(t),'color',[0.25 0.25 0.25,0.25],'linestyle','-','linewidth',1); 
hold on; plot(x,y_test(t),'color',[0.612 0.25 0.45],'linestyle','-','linewidth',1);
ylabel('Firing Rate'); 
yyaxis right
hold on; plot(x,x_test(t),'color',[0 0.25 0.612 0.25],'linestyle','-','linewidth',1);
ylabel('DFF','Interpreter','latex'); 
set(gca,'YColor','k');
xlabel('Time (s)');
title('Exampled Data from testing Set');

%testing data other methods
figure('units','centimeters','position',[13 13 15 4]); 
hold on; plot(x,spikes_probe_test(t),'color',[0.25 0.25 0.25,0.25],'linestyle','-','linewidth',1); 
plot(x,y_glm(t),'color',[0.612 0.25 0.45],'linestyle','-','linewidth',1); 
ylabel('Firing Rate'); 
yyaxis right
plot(x,y_glm_lr(t),'color',[0.25 0.612 0.45],'linestyle','-','linewidth',1);
plot(x,y_lr(t),'color',[0.1 0.2 0.45],'linestyle','-','linewidth',1);
ylabel('Firing Rate'); 
hold on; plot(x,y_lr(t),'color',[0 0.25 0.612 0.25],'linestyle','-','linewidth',1);
set(gca,'YColor','k');
xlabel('Time (s)');
legend({'fr','glm','lr_glm','lr'},'Interpreter','none');
title('Alternative Methods');

%autocorrelation and cross correlation with true
figure; hold on; 
n = numel(y_test);
data = cat(2,spikes_probe_test(1:n),trace_probe_test(1:n),y_test,y_glm(1:n),y_glm_lr(1:n),y_lr(1:n));
data = data-nanmean(data,1);
label = {'fr','dff','narx','glm','lr_glm','lr'};
t = [-60:60]/30;
for i = 1:numel(label)
   subplot(3,2,i);    hold on;
   plot(t,xcorr(data(:,1),data(:,i),(numel(t)-1)/2,'normalized'),'linewidth',2,'color','k');   
   ylabel('rho')
   yyaxis right
   plot(t,xcorr(data(:,i),data(:,i),(numel(t)-1)/2,'normalized'),'linewidth',2,'color',[0.25 0.612 0.45]);
   plot([0,0],[-0.2 1],'color','r','linestyle',':','linewidth',2);   
   set(gca,'YColor','k');
   ylabel('rho')
   title(label{i},'Interpreter','none');
%    legend('xcorr','auto');
   xlabel('lag (s)')
end
sgtitle('Xcorr (black) and Autocorr (green) Across Methods');

handles = get(groot, 'Children');
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\SupplementalFigures\Deconvolution\';
for i =1:numel(handles)
    saveas(handles(i),[savedir,sprintf('deconv_examples_probe2_rec4_500uM_%d.fig',i)])
end
close all

%% write the csd code
%write general code for figures tomorrow

%% Compare DFF to Neural across layers

%% implement a feedforward model structure
% loop through different parameters and identify best network on the
% training data (closed loop). 



%%
%overall figures to plot across all within recordings
%autocorrelation of imaging vs deconovlved vs ephys
%average kernel 
%correlation down 
%% Across all recordings, sweep layers (Figure 2)
%load the csd


%sweep layers (and their combinations) using each method



%plot the results 




%% Choose a the optimal average depth across all probes and run the rest of the comparisons. (Figure 3)




params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 500]; %depth from surface of probe
params.radius = 1; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
params.splittype = 'half'; %type of training/test split to use (current just first and second half of recording
[dff,~,~] = CompileData_deconvolution(dff_list(4),spike_opts_list(4),params);
%loop through depths 
depths = [0:250:900; 250:250:1000]';
peakrho = [];
for cur_depth = 1:size(depths,1)
    params.depth = depths(cur_depth,:);
    [~,st,~] = CompileData_deconvolution([],spike_opts_list(4),params);
    %loop through probes and best xcorr (wihin a 1 sec)
    for cur_probe = 1:4
        trace_probe = (dff{1}(:,cur_probe));
        spikes_probe = (st{1}(:,cur_probe));
        [a,~] = xcorr(spikes_probe-nanmean(spikes_probe),trace_probe-nanmean(trace_probe),120,'normalized');
        peakrho(cur_depth,cur_probe) = max(a);
    end
end
open peakrho

%%














%% COMPARE deconvolution methods across depths
%to do: load csd and compute true depth with that 
%make depths special per probe: verticalDepth(depth_subset(mode(source)),spike_opts.surface_offset,csd_opts.probe_angle(cur_probe));
% compile data
depths = [0:250:900; 250:250:1000]';
type_depth = cell(1,size(depths,1));
stats_depth = cell(1,size(depths,1));
for cur_depth = 1:size(depths,1)
    params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
    params.depth = depths(cur_depth,:); %depth from surface of probe
    params.radius = 2; %pixel radius around probe tip    
    params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
    params.splittype = 'half'; %type of training/test split to use (current just first and second half of recording
    if cur_depth==1
        [dff,st,~] = CompileData_deconvolution(dff_list(1:2),spike_opts_list(1:2),params);
    else %don't reload imaging
        [~,st,~] = CompileData_deconvolution([],spike_opts_list(2),params);
    end

    %train on first half of recordings
    dff_train = cellfun(@(x) x(1:floor(size(x,1)/2),:),dff,'UniformOutput',0);
    st_train = cellfun(@(x) x(1:floor(size(x,1)/2),:),st,'UniformOutput',0);
    dff_test = cellfun(@(x) x(floor(size(x,1)/2)+1:end,:),dff,'UniformOutput',0);
    st_test = cellfun(@(x) x(floor(size(x,1)/2)+1:end,:),st,'UniformOutput',0);
    trained_opts = Deconvolve_Train(dff_train,st_train,'all');

    % Test within each probe per recroding
    probe_idx = repmat(1:4,numel(dff_test),1)';
    rec_idx = repmat(1:numel(dff),4,1);
    train_idx = [rec_idx(:),probe_idx(:)];
    %test on the same data (second have of the recording)
    test_idx = train_idx; 
    stats = cell(size(train_idx,1),4);
    for i = 1:size(train_idx,1)    
        fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
        stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'narx',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));   
    end %rec
    %concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
    type = repmat({'narx','lr_gcamp','lr_glm','glm'},size(stats,1),1);
    type_depth{cur_depth} = cat(1,type(:));
    stats_depth{cur_depth} = cat(1,stats{:});
end

%% Choose Depth
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = 500; %depth from surface of probe
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
params.splittype = 'half'; %type of training/test split to use (current just first and second half of recording
[dff,st,fig_handles] = CompileData_deconvolution(dff_list(1:2),spike_opts_list(1:2),params);

%train on first half of recordings
dff_train = cellfun(@(x) x(1:floor(size(x,1)/2),:),dff,'UniformOutput',0);
st_train = cellfun(@(x) x(1:floor(size(x,1)/2),:),st,'UniformOutput',0);
dff_test = cellfun(@(x) x(floor(size(x,1)/2)+1:end,:),dff,'UniformOutput',0);
st_test = cellfun(@(x) x(floor(size(x,1)/2)+1:end,:),st,'UniformOutput',0);
trained_opts = Deconvolve_Train(dff_train,st_train,'all');

% Test within each probe per recroding
probe_idx = repmat(1:4,numel(dff_test),1)';
rec_idx = repmat(1:numel(dff),4,1);
train_idx = [rec_idx(:),probe_idx(:)];
%test on the same data (second have of the recording)
test_idx = train_idx; 
stats = cell(size(train_idx,1),4);
for i = 1:size(train_idx,1)    
    frpintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
    stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'narx',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));   
end %rec
%concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
type = repmat({'narx','lr_gcamp','lr_glm','glm'},size(stats,1),1);
type_within = cat(1,type(:));
stats_within = cat(1,stats{:});


%% Test across sites within recordings
%each train will be used 3 times (other locations)
temp = repmat(1:4,3,1);
probe_idx = repmat(temp(:)',numel(dff_test),1)';
rec_idx = repmat(1:numel(dff),numel(temp),1);
train_idx = [rec_idx(:),probe_idx(:)];

temp = repmat(1:4,4,1)';
temp(logical(eye(size(temp)))) = []; %remove auto comparisons
probe_idx = repmat(temp(:)',numel(dff_test),1)';
rec_idx = repmat(1:numel(dff),numel(temp),1);
test_idx = [rec_idx(:),probe_idx(:)];

stats = cell(size(train_idx,1),4);
for i = 1:size(train_idx,1)    
    frpintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
    stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'narx',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));   
end %rec
%concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
type = repmat({'narx','lr_gcamp','lr_glm','glm'},size(stats,1),1);
type_acrosssites = cat(1,type(:));
stats_acrosssites = cat(1,stats{:});

%% Test across days, same site, same animal
%get recording ids 
mousenum = MouseNumFromPath(dff_list,'Mouse_');
unique_mice = unique(mousenum);

%each train will be used onces
probe_idx = repmat(1:4,numel(dff_test),1)';
rec_idx = repmat(1:numel(dff),4,1);
train_idx = [rec_idx(:),probe_idx(:)];

%flip test days so diff than train days per mouse
test_idx = train_idx; 
for i = 1:numel(unique_mice)
    idx = find(mousenum==unique_mice(i));
    %flip the rec numbers
    for j = 1:numel(idx)
        test_idx(train_idx(:,1)==idx(j),1)=setdiff(idx,j);
    end
end

stats = cell(size(train_idx,1),4);
for i = 1:size(train_idx,1)    
    frpintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
    stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'narx',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
    stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));   
end %rec
%concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
type = repmat({'narx','lr_gcamp','lr_glm','glm'},size(stats,1),1);
type_acrossrec = cat(1,type(:));
stats_acrossrec = cat(1,stats{:});

%% Make figures


%to do: 
%choose on recording and make detailed figures and hone the network params
%Example deconvolution using all four methods for each probe
%Plotted using top 500uM;

%compare deconvolution quality by number of neurons (1:100); 

%sweep the depths in 100uM increments for the 'within fits'
%plot the quality of fit within rec/site as a function of depth per
%deconvolution type

%For 'best depth':
%compare between sites and between recs
%exmample traces of deconvolution methods on 1 recording
%example traces of same test location with different train location in same
%recording
%example traces of same test location with different train day

%basic statistics to plot
% plot the csd
% plot histogram of unit depths per probe
% plot the number of neurons. 
% add ignore bad ones ot the compile methdos


% 
% 
% %% make plots comparing the methods
% %compare correlation to withheld data
% figure; hold on; 
% data = arrayfun(@(n) diag(stats{n}.rho), 1:numel(stats),'UniformOutput',0);
% data = [data{:}];
% bar(nanmean(data,1))
% 
% %compare generalizability across areas
% figure; hold on; 
% idx = eye(size(stats{1}.err));
% data = arrayfun(@(n) stats{n}.err(~idx), 1:numel(stats),'UniformOutput',0);
% data = [data{:}];
% bar(nanmean(data,1))



