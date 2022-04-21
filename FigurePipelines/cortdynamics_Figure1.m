%Figure 1 
%Camden MacDowell

%% Probe tracing and visualization
fp = fig_params_cortdynamics;
%uncomment to run (these can take a bit and are older functions)
% AnatomicalCorrelellogram(EphysPath,probe_ccf_path,allen_atlas_path,mua_flag,fp,cur_probe)
% PlotCombinedProbeTrajectories(histo_wrkspace_path,spike_opts_path)

% Get a summary of the in/out correlations across recordings per area
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ValidationProbeLocation';
Plot_AnatomicalCorrelationInOut(1,'general');
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'In_OutCorrelations',savedir,0); close all



savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ValidationProbeLocation';
for i = 1:6
    PlotExampleCorrellelogram(i)
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
    set(gca,'CLim',[0 0.025])    
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ActivityCorrelelogram',savedir,0); close all


%% Plot example raw data 
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ExampleData';
if ~exist(savedir,'dir'); mkdir(savedir); end

cur_rec = 1;
tp = [1200 2100];
PlotExampleEphys(1,1,tp)
fp.FigureSizing(gcf,[3 2 16 6],[2 10 30 15])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleEphys',savedir,0); close all

PlotExampleImagingTrace(cur_rec,tp)
fp.FigureSizing(gcf,[3 2 16 3],[2 10 40 15]); box on
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleImaging',savedir,0); close all

%plot motif onsets
motif_num = [3,5,14];
PlotExampleMotifWeightings(cur_rec,motif_num,tp); 
fp.FigureSizing(gcf,[3 2 16 1],[2 10 40 15])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleMotifOnsets',savedir,0); close all

%plot motifs
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Motifs';
if ~exist(savedir,'dir'); mkdir(savedir); end
PlotBasisMotifs(savedir) 

%plot the average response of ephys areas to motifs
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\EphysResponse';
if ~exist(savedir,'dir'); mkdir(savedir); end
Plot_MotifEphysActivity(1,1:14,1,'mean',savedir)

%% Load subspace data
data = LoadSubspaceData('in_grouped');
dataout = LoadSubspaceData('out_grouped');

%% ridge regression figures
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\RidgeRegressionFigures';
if ~exist(savedir,'dir'); mkdir(savedir); end
%for IN models
RidgeRegressionExampleFigures(data,'VIS',[5,3])
RidgeRegressionCombinedFigures(data)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'RegressionExample_IN',savedir,0); close all

%for OUT models
RidgeRegressionExampleFigures(dataout,'VIS',[5,3])
RidgeRegressionCombinedFigures(dataout)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'RegressionExample_OUT',savedir,0); close all



%next to add: 
%plot showing the correlation between data in and data out... explain in
%text that this is not required and highlight it's difference in strength
%between different areas. 
%plot showing the decodability in average |trial-to-trial activity| across all neurons between each pair of motifs.  
%the permutation tests
%by end of next week | knock out figure 2 and have explored 3


%%




















