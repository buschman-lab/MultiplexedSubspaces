function GutcheckDeconvolution(stack,gp)
%stack is a masked, NxPixel DFF trace
if nargin <2; gp = general_params; end
fp = fig_params;
gp.w_deconvolution = 'only_filter'; %original
[~, nanpxs_orig, data_orig, ~] = ProcessAndSplitData(stack,[],gp);

gp.w_deconvolution = 'filter_thresh'; %thresholded
[~, nanpxs_thresh, data_thresh, data_thresh_test] = ProcessAndSplitData(stack,[],gp);

gp.w_deconvolution = 'lucric'; %deconvolved with lucy richardson
[~, nanpxs_lucrid, data_lucric, data_lucric_test] = ProcessAndSplitData(stack,[],gp);

%Compare the dynamics just using the first chunk
figure('position',[680    79   677   899]); hold on; 
sgtitle('Comparing Methods: 2 min chunk')
subplot(3,1,1); imagesc(squeeze(data_orig(:,:,1))); set(gca,'xtick',[],'ytick',[]); ylabel('orig'); colorbar; fp.FormatAxes(gca);
subplot(3,1,2); imagesc(squeeze(data_thresh(:,:,1))); set(gca,'xtick',[],'ytick',[]); ylabel('threshold (1STD)'); colorbar; fp.FormatAxes(gca);
subplot(3,1,3); imagesc(squeeze(data_lucric(:,:,1))); set(gca,'xticklabel',round(get(gca,'xtick')/gp.fps,0),'ytick',[]); ylabel('lucy-goosey'); colorbar; xlabel('time (s)'); fp.FormatAxes(gca);
colormap magma;

%zoom in on the the first 10 seconds
figure('position',[680    79   677   899]); hold on; 
idx =(1:round(gp.fps*30,0));
sgtitle('Comparing Methods: Zoomed In')
subplot(3,1,1); imagesc(squeeze(data_orig(:,idx,1))); set(gca,'xtick',[],'ytick',[]); ylabel('orig'); colorbar; fp.FormatAxes(gca);
subplot(3,1,2); imagesc(squeeze(data_thresh(:,idx,1))); set(gca,'xtick',[],'ytick',[]); ylabel('threshold (1STD)'); colorbar; fp.FormatAxes(gca);
subplot(3,1,3); imagesc(squeeze(data_lucric(:,idx,1))); set(gca,'xticklabel',round(get(gca,'xtick')/gp.fps,0),'ytick',[]); ylabel('lucy-goosey'); colorbar; xlabel('time (s)'); fp.FormatAxes(gca);
colormap magma;

%Compare the discovered motifs on thresh vs deconvolved
fprintf('\nWorking on Thresholded Data\n');
[W_thresh, H_thresh, stats_train_thresh, stats_test_thresh] = FitandValidateMotifs(squeeze(data_thresh(:,:,1)),squeeze(data_thresh_test(:,:,1)),gp,0);    
fprintf('\nWorking on Deconvolved Data\n');
[W_lucric, H_lucric, stats_train_lucric, stats_test_lucric] = FitandValidateMotifs(squeeze(data_lucric(:,:,1)),squeeze(data_lucric_test(:,:,1)),gp,0);

%Show a figure plotting all of their flattened Ws, along with the PEV and other stats. 

%organize both groups by the motif with the highest expvar

%then I guess keep ploting until you run out of one and then continue with
%the rest


end