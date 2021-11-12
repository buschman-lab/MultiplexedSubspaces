function ImpactOfFluorescenceNormalization(mouse)
%Camden MacDowell - timeless
%grab the data for one animal/one recording
% mouse = {'331_1','331_2','332_1','332_2','334_1','334_2'}
fp = fig_params_deconvolutionpaper;
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Normalization\';
switch mouse
    case '332_1'
        baseline_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\BaselineRestingState\Mouse332_06_02_2021\Mouse332_RestingState_06_02_2021_1'}); %baseline raw stack
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Mouse332_RestingState_NP_06_07_2021_1'}); %cranio raw stack
    case '332_2'
        baseline_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\BaselineRestingState\Mouse332_06_02_2021\Mouse332_RestingState_06_02_2021_1'}); %baseline raw stack
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_08_2021\Mouse332_RestingState_NP_06_08_2021_1'}); %cranio raw stac
    case '331_1'
        baseline_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\BaselineRestingState\Mouse331_06_04_2021\Mouse331_RestingState_06_04_2021_1'}); %baseline raw stack
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_11_2021\Mouse331_RestingState_NP_06_11_2021_1'}); %cranio raw stac        
    case '331_2'
        baseline_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\BaselineRestingState\Mouse331_06_04_2021\Mouse331_RestingState_06_04_2021_1'}); %baseline raw stack
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_12_2021\Mouse331_RestingState_NP_06_12_2021_1'}); %cranio raw stac        
    case '334_1'
        baseline_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\BaselineRestingState\Mouse334_06_03_2021\Mouse334_RestingState_06_03_2021_1'}); %baseline raw stack
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_09_2021\Mouse334_RestingState_NP_06_09_2021_1'}); %cranio raw stac  
    case '334_2'
        baseline_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\BaselineRestingState\Mouse334_06_03_2021\Mouse334_RestingState_06_03_2021_1'}); %baseline raw stack
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_10_2021\Mouse334_RestingState_NP_06_10_2021_1'}); %cranio raw stac 
end

processed_path = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\'; %for probe coords
[~,fn_header] = fileparts(craniotomy_file{1});
probe_coords = load([processed_path,erase(fn_header,'_MMStack_Pos0.ome_stack'),'dff_combined.mat'],'probe_coords'); 
probe_coords=probe_coords.probe_coords;
fig_header = erase(fn_header,'_MMStack_Pos0.ome_stack');

% Plot the variance of the baseline recording (processed as dff)
dff = GrabAFewMinutes(baseline_file,1); %dff_flag==1-->hemo 2=nohemo
plotVarianceImage(dff,'std'); colorbar
saveCurFigs(gcf,'-png',[fig_header,'baseline_dff_std'],savedir,1); close all

% Plot the variance of the recording with craniotomies (dff and dfs)
[stack,opts] = GrabAFewMinutes(craniotomy_file,0); %dff_flag==1-->hemo 2=nohemo 0=nothing
opts.method = 'movingavg';
dff = makeDFF(stack, opts, 'dff', opts.method_window);
opts.method = 'zscore';
dfs = makeDFF(stack, opts, 'dff', opts.method_window);
plotVarianceImage(dff,'std'); colorbar; %caxis([0 0.25])
saveCurFigs(gcf,'-png',[fig_header,'cranio_dff_std'],savedir,1); close all
plotVarianceImage(dfs,'std'); colorbar; caxis([0 1])
saveCurFigs(gcf,'-png',[fig_header,'cranio_dfs_std'],savedir,1); close all

% Plot the correlation between hemispheres (use dff)
[rho_in, rho_out] = CraniotomyContraCorrelation(dff,probe_coords);
data = [rho_out,cat(1,rho_in,NaN(numel(rho_out)-numel(rho_in),1))];
figure; CompareViolins(data',fp,'sidebyside',1,'col',{[0.2 0.2 0.2],[0.91 0.41 0.17]});
ylabel('Rho_z','Interpreter','tex');
set(gca,'units','centimeters','position',[3 2 2.5 4],'xlim',[0.75 1.25]); fp.FormatAxes(gca); 
[~,p] = kstest2(rho_in,rho_out); title(sprintf('KS p=%0.2g',p),'fontweight','normal');
saveCurFigs(gcf,'-svg',[fig_header,'correlation_between_hemis'],savedir,1); close all

% plot the histograms comparing pixels values in each hemisphere
PlotHemisphereVariance(dff,'all',100,[-2 2]); 
set(gca,'units','centimeters','position',[2 2 4 4]); fp.FormatAxes(gca); 
saveCurFigs(gcf,'-svg',[fig_header,'variance_between_hemis_dff'],savedir,1); close all
PlotHemisphereVariance(dfs,'all',100,[-3 3]); 
set(gca,'units','centimeters','position',[2 2 4 4]); fp.FormatAxes(gca); 
saveCurFigs(gcf,'-svg',[fig_header,'variance_between_hemis_dfs'],savedir,1); close all

% plot the violins comparing each craniotomy std 
[activity_in,activity_contra]=CompareROIVariance(dff,probe_coords,'pxlstd');
figure; CompareViolins(activity_in',fp,'col',{[0.91 0.41 0.17],[0.91 0.41 0.17],[0.91 0.41 0.17],[0.91 0.41 0.17]});
set(gca,'units','centimeters','position',[2 2 6 4]); fp.FormatAxes(gca); ylabel('std') 
saveCurFigs(gcf,'-svg',[fig_header,'across_cranio_var_dff'],savedir,1); close all

figure; CompareViolins(activity_contra',fp);
set(gca,'units','centimeters','position',[2 2 6 4]); fp.FormatAxes(gca); ylabel('std') 
saveCurFigs(gcf,'-svg',[fig_header,'across_cranio_var_dff'],savedir,1); close all

% plot the violins comparing each craniotomy all 
[activity_in,activity_contra]=CompareROIVariance(dff,probe_coords,'all');
figure; CompareViolins(activity_in',fp,'col',{[0.91 0.41 0.17],[0.91 0.41 0.17],[0.91 0.41 0.17],[0.91 0.41 0.17]});
set(gca,'units','centimeters','position',[2 2 6 4]); fp.FormatAxes(gca); ylabel('dff') 
saveCurFigs(gcf,'-svg',[fig_header,'across_cranio_all_dff'],savedir,1); close all

figure; CompareViolins(activity_contra',fp);
set(gca,'units','centimeters','position',[2 2 6 4]); fp.FormatAxes(gca); ylabel('dff') 
saveCurFigs(gcf,'-svg',[fig_header,'across_cranio_all_dff'],savedir,1); close all

% plot the violins comparing each craniotomy all 
[activity_in,activity_contra]=CompareROIVariance(dfs,probe_coords,'all');
figure; CompareViolins(activity_in',fp,'col',{[0.91 0.41 0.17],[0.91 0.41 0.17],[0.91 0.41 0.17],[0.91 0.41 0.17]});
set(gca,'units','centimeters','position',[2 2 6 4]); fp.FormatAxes(gca); ylabel('dff') 
saveCurFigs(gcf,'-svg',[fig_header,'across_cranio_all_dfs'],savedir,1); close all

figure; CompareViolins(activity_contra',fp);
set(gca,'units','centimeters','position',[2 2 6 4]); fp.FormatAxes(gca); ylabel('dfs') 
saveCurFigs(gcf,'-svg',[fig_header,'across_cranio_all_dfs'],savedir,1); close all

% plot the variance in the craniotomies vs the variance in the other hemi
[activity_in,activity_contra]=CompareROIVariance(dff,probe_coords,'all');
PlotCraniotomyVariance(activity_in,activity_contra,100,[-1 1]);
set(gca,'units','centimeters','position',[2 2 3 4]); fp.FormatAxes(gca); ylabel('%values') 
saveCurFigs(gcf,'-svg',[fig_header,'cranio_contra_variance_dff'],savedir,1); close all

%get the pairwise F statistic between cranio sites and between contra
idx = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
fstat_in_dff = NaN(size(idx,1),1);
fstat_out_dff = NaN(size(idx,1),1);
for i = 1:size(idx,1)
    [~,~,~,stats] = vartest2(activity_in(:,idx(i,1)),activity_in(:,idx(i,2)));
    fstat_in_dff(i) = stats.fstat;
    [~,~,~,stats] = vartest2(activity_contra(:,idx(i,1)),activity_contra(:,idx(i,2)));
    fstat_out_dff(i) = stats.fstat;
end
%get the between as well
[~,~,~,stats] = vartest2(activity_contra(:),activity_in(:));
fstat_between_dff = stats.fstat;

% plot the variance in the craniotomies vs the variance in the other hemi
[activity_in,activity_contra]=CompareROIVariance(dfs,probe_coords,'all');
PlotCraniotomyVariance(activity_in,activity_contra,100,[-3 3]);
set(gca,'units','centimeters','position',[2 2 4 4]); fp.FormatAxes(gca); ylabel('%values') 
saveCurFigs(gcf,'-svg',[fig_header,'cranio_contra_variance_dfs'],savedir,1); close all

%get the pairwise F statistic between cranio sites and between contra
idx = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
fstat_in_dfs = NaN(size(idx,1),1);
fstat_out_dfs = NaN(size(idx,1),1);
for i = 1:size(idx,1)
    [~,~,~,stats] = vartest2(activity_in(:,idx(i,1)),activity_in(:,idx(i,2)));
    fstat_in_dfs(i) = stats.fstat;
    [~,~,~,stats] = vartest2(activity_contra(:,idx(i,1)),activity_contra(:,idx(i,2)));
    fstat_out_dfs(i) = stats.fstat;
end
%get the between as well
[~,~,~,stats] = vartest2(activity_contra(:),activity_in(:));
fstat_between_dfs = stats.fstat;

save([savedir fig_header, '_stats.mat'],'fstat_in_dfs','fstat_out_dfs','fstat_in_dff','fstat_out_dff','fstat_between_dff','fstat_between_dfs','rho_in','rho_out');


end

function [stack,opts] = GrabAFewMinutes(file_list,dff_flag)
    %grab multiple preprocessed stacks (but not the full hour to save time)
    for i = 1:15 %five should be enough its  a little over 20 min at 30Hz
       [file_path,fn] = fileparts(file_list{1}); 
       if i==1
           temp_stack = load(file_list{1},'stack','opts');
           opts = temp_stack.opts;
           stack = temp_stack.stack;
       else
           temp_stack = load([file_path filesep erase(fn,'.ome_stack') sprintf('_%d.ome_stack.mat',i-1)],'stack');
           stack = cat(3,stack,temp_stack.stack);
       end
    end
    
    if dff_flag==1 %corrected
       stack = HemodynamicCorrection(stack, opts);
    elseif dff_flag==2 %uncorrected
       [~,stack] = HemodynamicCorrection(stack, opts);
    else %filter to remove
        stack = filterstack(stack, 30, [0.1 4], 'lowpass', 1, 0);
    end
    
    mask = repmat(imresize(opts.mask,[68 68]),1,1,size(stack,3)); 
    stack(mask==0)=0; %correct in case this is yet to be included in the preprocessing
end




%Tonight/Tomorrow morning: 
%retrace probe locations (do them perfectly)
%make the figures for each animal
%then write the function for all of the ephys as the final plot
%write function to plot on other peoples available data 


% %Compare Ephys 
% load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\restingstate_processed_fn.mat','spike_opts_list') 
% params.mua = 0; %1= use both 'good' and 'mua' units. 0 = just 'good'
% params.depth = [0 600]; %depth from surface of probe
% params.radius = 2; %pixel radius around probe tip    
% params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
% params.bindata = 0; %1=temporally bin data (i.e. to frame rate 15 from 30)
% [~,st,n_neurons] = CompileData_deconvolution([],spike_opts_list(3),params); 
% 
% %convert to firing rate (hZ) and normalize for the number of neurons and get the std
% fr_std = arrayfun(@(n) nanstd(30*(st{n}./n_neurons(n,:)),[],1),1:numel(st),'UniformOutput',0);
% fr_std = cat(1,fr_std{:});
% 
% grp = [1:4].*ones(numel(st),1);
% anovan(fr_std(:),grp(:));
% 




