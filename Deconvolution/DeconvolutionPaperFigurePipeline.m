%% Deconvolution Figure Pipeline
%Camden MacDowell - timeless

%merge files and resave
%include plot of full trained vs generalized ff
%decide how to plot the firing rate distribution for each type averaged across all
%make into figures
%get the drift of each neuron for each recording

%% Example Traces on trained or withheld data
% savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\ExampleTracesTraining';
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\ExampleTracesWithheld';
fp = fig_params_deconvolutionpaper; 
col = [fp.c_none; fp.c_lr; fp.c_glm ;fp.c_ff];
% fn = GrabFiles('deconvolved_traces_trained_datastd.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'}); %trained
% fn = GrabFiles('deconvolved_traces_trained_data_10pvalidationAndGLMinterceptstd.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'}); savefn = 'trained'; %trained
fn = GrabFiles('deconvolved_traces_withheld_datastd.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'}); savefn = 'withheld'; %withheld
data = load(fn{1});
labels = {'None','LR','GLM','fNN'};
%loop through each rec and each probe and generate the figure
n_rec = size(data(1).st_test,2);
n_probe = size(data(1).st_test{1},2);
for cur_rec = 1:n_rec
    for cur_probe = 1:n_probe
       %concatenate data for this probe   
       close all
       ypred = [data(1).ypred_none{cur_rec}(:,cur_probe),data(1).ypred_lr{cur_rec}(:,cur_probe),...
           data(1).ypred_glm{cur_rec}(:,cur_probe),data(1).ypred_ff{cur_rec}(:,cur_probe)];
       y = data(1).st_test{cur_rec}(:,cur_probe);
       x = (1:size(ypred,1))/30;
       subdur=20; %duration of the subplot in seconds
       posdur=190; %position of the subplot in seconds
       %plot each individually, including true, for 5 minutes 
       for cur_method = 1:size(ypred,2)+1        
           figure; hold on; 
           %add zoom box
           if cur_method==size(ypred,2)+1
              yvals = [floor(min(y(1:600*30))),ceil(max(y(1:600*30)))];
           else
              yvals = [floor(min(ypred(1:600*30,cur_method))),ceil(max(ypred(1:600*30,cur_method)))];
           end           
           r=rectangle('position',[posdur,yvals(1),subdur,sum(abs(yvals))],'FaceColor',[0.25 0.25 0.25 0.25],'EdgeColor','none');
           plot([0 x(end)],[0 0],'color','k','linewidth',fp.line_width)
           if cur_method==size(ypred,2)+1
               p=plot(x,y,'linewidth',fp.line_width,'color','k');
           elseif cur_method >1 && cur_method <= size(ypred,2)
%                plot(x,y,'linewidth',fp.line_width,'color',[0.2 0.2 0.2, 0.1],'linewidth',1);
               p=plot(x,ypred(:,cur_method),'linewidth',fp.line_width,'color',col(cur_method,:));
           else
               p=plot(x,ypred(:,cur_method),'linewidth',fp.line_width,'color',col(cur_method,:));
           end
           if cur_method==1
               ylabel('\DeltaF/F'); 
           elseif cur_method == size(ypred,2)+1
               ylabel({'Normalized','firing rate'});
           else
               ylabel({'Predicted','firing rate'});
           end
           xlim([0 300]);
           xlabel('Time (s)')                                     
           fp.FormatAxes(gca); grid on
           if cur_method==size(ypred,2)+1
               fp.SetTitle(gca,'Original')
           else
               fp.SetTitle(gca,labels{cur_method})
           end           
           fp.FigureSizing(gcf,[3 2 20 3],[10 10 25 7])             
           %save off
           saveCurFigs(gcf,{'-dpng','-dsvg'},sprintf('Rec%d_P%d_%d',cur_rec,cur_probe,cur_method),savedir,0); 
           
           %zoom in, change size, add box
           delete(r); 
           xlim([posdur,posdur+subdur])
           fp.FigureSizing(gcf,[3 2 10 3],[10 10 15 7]) 
           p.LineWidth=1; box on 
           saveCurFigs(gcf,{'-dpng','-dsvg'},sprintf('Rec%d_P%d_%dzoom',cur_rec,cur_probe,cur_method),savedir,0); close
       end
    end
end

%% Within site/animals training data stats
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\training';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit to 600uM 
fn = GrabFiles('trainingstd.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1},{'lr_gcamp','glm','feedforward','none'});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,1,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%zoom in on glm and fNN plots for rho, err, and skew
sz = [2 2 2 4];
data_temp = cellfun(@(x) x(:,3:4),data,'UniformOutput',0);
fh_stats = makeViolinPlots(data_temp,1,fp,sz,col(3:4),labels(3:4));
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('statszoom'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all


%% Within site/animals
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\withinsites_animal';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit to 600uM 
fn = GrabFiles('within_compare\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{end-1},{'lr_gcamp','glm','feedforward','none'});%get the desired depth 
labels = {'None','LR','GLM','fNN'};
%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,1,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all

% figure comparing the distibutions shape
% figure comparing methods across depths

%% Across site within animal
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\xsite_within_animal';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];

%load fit Across Sites uM 
fn = GrabFiles('within_xsite\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1},{'lr_gcamp','glm','feedforward','none'});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,0,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all

%% Across recordings in same animal. same site
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\xrec_same_siteanimal';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit Across Sites uM 
fn = GrabFiles('xrec_samesitesstd\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,1,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all
%plot the statistical comparisons

sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all

%% Across animals. same site
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\xanimals_samesite';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit Across Sites uM 
fn = GrabFiles('xanimals\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,0,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all

%% train Across sites eval on testing data
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\train_xsite';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit Across Sites uM 
fn = GrabFiles('train_xsitesstd\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1},{'lr_gcamp','glm','feedforward','none'});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,1,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all


%% train Across sites eval on training data
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\train_xsite_evaltrain';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit Across Sites uM 
fn = GrabFiles('train_xsites_evaluatedontrainingdatastd\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1},{'lr_gcamp','glm','feedforward','none'});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,1,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all

%% train everything eval testing data
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\train_all';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit Across Sites uM 
fn = GrabFiles('train_xeverythingstd\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1},{'lr_gcamp','glm','feedforward','none'});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,1,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all

%% train everything evaluate training data
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\train_all_evaltrain';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit Across Sites uM 
fn = GrabFiles('train_xsites\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{1},{'lr_gcamp','glm','feedforward','none'});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%reverse skew
data{3} = -1*data{3};

%make plots for rho, err, and skew
sz = [2 2 6 8];
fh_stats = makeViolinPlots(data,1,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%plot the statistical comparisons
sz = [2 2 1.75 1.75];
handles = makeStatsPlots(data,fp,sz,labels);  
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('pvals'),savedir,0); close all

%plot small plots of the xcross correlation
sz = {[2 2 8 8],[2 2 4 4]};
fh_xcorr = makeXcorrPlots(xcorr_trace,fp,sz,col,labels);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('xcorr'),savedir,0); close all

%% plot the training across depths ON TESTING DATA
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\across_depth_test';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load data
fn = GrabFiles('within_compare\w*mean\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
n=7;
fn = fn(1:7);
data = cell(numel(fn),3);
for i = 1:numel(fn)
    data(i,:) = LoadResults(fn{i},{'lr_gcamp','glm','feedforward','none'});%get the desired depth     
end

labels = {'None','LR','GLM','fNN'};
%load depths
depths = load(fn{1},'depths');depths = nanmedian(depths.depths(1:n,:),2);
n_neurons = cellfun(@(x) load(x,'n_neurons'),fn,'UniformOutput',0);
n_neurons_avg = cellfun(@(x) nanmean(x.n_neurons(:)),n_neurons,'UniformOutput',1);
n_neurons_sem = cellfun(@(x) sem(x.n_neurons(:)),n_neurons,'UniformOutput',1);

statname = {'rho','err','s'};
for i = 1:3 %loop through statistics
    figure; hold on; 
    %equate their sizes
    maxlen = max(cellfun(@(x) size(x,1),data(:,i),'UniformOutput',1));
    temp = cellfun(@(x) cat(1,x, NaN(maxlen-size(x,1),size(x,2))),data(:,i),'UniformOutput',0);
    %loop through type
    for j=1:size(temp{1},2)
        data_method = cellfun(@(x) x(:,j), temp,'UniformOutput',0);
        data_method = [data_method{:}];

        y = nanmean(data_method,1)';
        x = depths;
        yerr = sem(data_method,1)';
        shadedErrorBar(x,y,yerr,'lineprops',{'linestyle','-','color',col{j}});    
        ylabel(statname{i})
        xlabel('Depth (uM)');              
    end
    fp.FormatAxes(gca); grid on
    fp.SetTitle(gca,statname{i})
    fp.FigureSizing(gcf,[2 4 8 6],[])      
    xlim([min(depths),max(depths)])    
    set(gca,'XTickLabelRotation',45)
end
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('depthcorr'),savedir,0); close all
figure; hold on; 
shadedErrorBar(depths,n_neurons_avg,n_neurons_sem,'lineprops',{'linestyle','-','color',[0.5 0.5 0.5]});    
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('depthneu'),savedir,0); close all

%% compare the superficial and deep
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\superficial_vs_deep';
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];

%load data
fn = GrabFiles('within_compare\w*mean\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn = fn(8:9);
type = {'lr_gcamp','glm','feedforward'};
data = cell(1,numel(type));
for j = 1:numel(type)        
   data{j} = LoadResultsByType(fn,type{j},'rho');
end

col = repmat([{fp.c_lr},{fp.c_glm},{fp.c_ff}],numel(fn),1);
col = col(:);    

%get the xvalues
offset=0.5; %gap between types
x = (repmat(1:numel(fn),numel(type),1)+[0:numel(fn)+offset:(numel(fn)+offset)*(numel(type)-1)]')';
x = x(:);
data = cat(1,data{:});

%plot
figure; hold on; 
vp = CompareViolins(data,fp,'plotspread',1,'xpos',x,'col',col);
xlim([0.5,x(end)+0.5])
ylabel('rho'); xlabel('Depth (um)')
set(gca,'xticklabels',{'0-600','600-1400'},'XTickLabelRotation',45)
pval = ranksum(data(1,:),data(2,:));
plot([x(1),x(2)],[0.7 0.7],'linewidth',1,'color','k'); 
text(nanmean(x(1:2)),0.7,sprintf('p=%0.2g',pval),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',12)
pval = ranksum(data(3,:),data(4,:));
plot([x(3),x(4)],[0.8 0.8],'linewidth',1,'color','k'); 
text(nanmean(x(3:4)),0.8,sprintf('p=%0.2g',pval),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',12)
pval = ranksum(data(5,:),data(6,:));
plot([x(5),x(6)],[0.9 0.9],'linewidth',1,'color','k'); 
text(nanmean(x(5:6)),0.9,sprintf('p=%0.2g',pval),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',12)
fp.FormatAxes(gca); grid on
fp.SetTitle(gca,'superficial vs deep')
fp.FigureSizing(gcf,[2 4 10 6],[])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('rho'),savedir,0); close all

% %% plot the training across depths ON training data
% savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\across_depth_TRAIN';
% fp = fig_params_deconvolutionpaper; 
% col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
% %load example traces and plot
% 
% %load data
% fn = GrabFiles('within_compare\w*mean\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
% n=7;
% fn = fn(1:7);
% data = cell(numel(fn),3);
% for i = 1:numel(fn)
%     data(i,:) = LoadResultsTrain(fn{i},{'lr_gcamp','glm','feedforward','none'});%get the desired depth     
% end
% 
% labels = {'None','LR','GLM','fNN'};
% %load depths
% depths = load(fn{1},'depths');depths = nanmedian(depths.depths(1:n,:),2);
% n_neurons = cellfun(@(x) load(x,'n_neurons'),fn,'UniformOutput',0);
% n_neurons_avg = cellfun(@(x) nanmean(x.n_neurons(:)),n_neurons,'UniformOutput',1);
% n_neurons_sem = cellfun(@(x) sem(x.n_neurons(:)),n_neurons,'UniformOutput',1);
% 
% statname = {'rho','err','s'};
% for i = 1:3 %loop through statistics
%     figure; hold on; 
%     %equate their sizes
%     maxlen = max(cellfun(@(x) size(x,1),data(:,i),'UniformOutput',1));
%     temp = cellfun(@(x) cat(1,x, NaN(maxlen-size(x,1),size(x,2))),data(:,i),'UniformOutput',0);
%     %loop through type
%     for j=1:size(temp{1},2)
%         data_method = cellfun(@(x) x(:,j), temp,'UniformOutput',0);
%         data_method = [data_method{:}];
% 
%         y = nanmean(data_method,1)';
%         x = depths;
%         yerr = sem(data_method,1)';
%         shadedErrorBar(x,y,yerr,'lineprops',{'linestyle','-','color',col{j}});    
%         ylabel(statname{i})
%         xlabel('Depth (uM)');              
%     end
%     fp.FormatAxes(gca); grid on
%     fp.SetTitle(gca,statname{i})
%     fp.FigureSizing(gcf,[2 4 8 6],[])      
%     xlim([min(depths),max(depths)])    
%     set(gca,'XTickLabelRotation',45)
% end
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('depthcorr'),savedir,0); close all
% figure; hold on; 
% shadedErrorBar(depths,n_neurons_avg,n_neurons_sem,'lineprops',{'linestyle','-','color',[0.5 0.5 0.5]});    
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('depthneu'),savedir,0); close all
% 
% %% compare the superficial and deep
% savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\superficial_vs_deep_train';
% fp = fig_params_deconvolutionpaper; 
% col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
% 
% %load data
% fn = GrabFiles('within_compare\w*mean\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
% fn = fn(8:9);
% type = {'lr_gcamp','glm','feedforward'};
% data = cell(1,numel(type));
% for j = 1:numel(type)        
%    data{j} = LoadResultsByType(fn,type{j},'rho');
% end
% 
% col = repmat([{fp.c_lr},{fp.c_glm},{fp.c_ff}],numel(fn),1);
% col = col(:);    
% 
% %get the xvalues
% offset=0.5; %gap between types
% x = (repmat(1:numel(fn),numel(type),1)+[0:numel(fn)+offset:(numel(fn)+offset)*(numel(type)-1)]')';
% x = x(:);
% data = cat(1,data{:});
% 
% %plot
% figure; hold on; 
% vp = CompareViolins(data,fp,'plotspread',1,'xpos',x,'col',col);
% xlim([0.5,x(end)+0.5])
% ylabel('rho'); xlabel('Depth (um)')
% set(gca,'xticklabels',{'0-600','600-1400'},'XTickLabelRotation',45)
% pval = ranksum(data(1,:),data(2,:));
% plot([x(1),x(2)],[0.7 0.7],'linewidth',1,'color','k'); 
% text(nanmean(x(1:2)),0.7,sprintf('p=%0.2g',pval),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',12)
% pval = ranksum(data(3,:),data(4,:));
% plot([x(3),x(4)],[0.8 0.8],'linewidth',1,'color','k'); 
% text(nanmean(x(3:4)),0.8,sprintf('p=%0.2g',pval),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',12)
% pval = ranksum(data(5,:),data(6,:));
% plot([x(5),x(6)],[0.9 0.9],'linewidth',1,'color','k'); 
% text(nanmean(x(5:6)),0.9,sprintf('p=%0.2g',pval),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',12)
% fp.FormatAxes(gca); grid on
% fp.SetTitle(gca,'superficial vs deep')
% fp.FigureSizing(gcf,[2 4 10 6],[])
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('rho'),savedir,0); close all

%% Compare within methods, across train/fit types
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\methods_acrosscomparisons';
fp = fig_params_deconvolutionpaper; 
fn = cell(9,1);
temp = GrabFiles('trainingstd.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn(1) = temp(2);
fn(2) = GrabFiles('within_comparedepthsstd8.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn(3) = GrabFiles('within_xsite\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn(4) = GrabFiles('xrec_samesites\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn(5) = GrabFiles('xanimals\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
temp = GrabFiles('train_xsites\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn(6) = temp(2);
fn(7) = GrabFiles('train_xsites_ev\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
temp = GrabFiles('train_xeverything\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
fn(8) = temp(2);
fn(9) = GrabFiles('train_xeverything_ev\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
type = {'lr_gcamp','glm','feedforward'};
titlestr = {'Correlation','MSE','Skew'};
ylabelstr = {'Rho-z','MSE','\Delta Skew'};
statname = {'rho','err','s'};
labels = {'LR','GLM','fNN'};
for i = 1:numel(statname)

    data = cell(1,numel(type));
    for j = 1:numel(type)        
       data{j} = LoadResultsByType(fn,type{j},statname{i});
    end

    col = repmat([{fp.c_lr},{fp.c_glm},{fp.c_ff}],numel(fn),1);
    col = col(:);    

    %get the xvalues
    offset=1; %gap between types
    x = (repmat(1:numel(fn),numel(type),1)+[0:numel(fn)+offset:(numel(fn)+offset)*(numel(type)-1)]')';
    x = x(:);
    data = cat(1,data{:});
    
    %reverse skew
    if i==3
        data = -1*data;
    end

    %plot
    figure; hold on; 
    vp = CompareViolins(data,fp,'plotspread',0,'xpos',x,'col',col);
    xlim([0.5,x(end)+0.5])
    ylabel(ylabelstr{i}); xlabel('Method')
    set(gca,'xtick',[2.5,7.5,12.5],'xticklabels',labels)
    fp.FormatAxes(gca); grid on
    fp.SetTitle(gca,titlestr{i})
    fp.FigureSizing(gcf,[2 2 12 6],[])
    
    if i==2
        yval = get(gca,'ylim');
        ylim([0 yval(2)]);
    elseif i==3
        xval = get(gca,'xlim');
        plot(xval,[0 0],'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--')
    end    
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('stats'),savedir,0); close all

%% gather data
function data = LoadResultsByType(fn,comptype,statname)
%each cell is a different method, each column is a different statistic
%input type allows either multiple method or multiple
    
temp = cell(1,numel(fn));
for i = 1:numel(fn)
    res = load(fn{i});
    if nargin <4
       deconv_type = {'lr_gcamp','glm','feedforward','none'}; 
    end
    deconv_type = repmat(deconv_type,size(res.test_idx,1),1);
    if ~isempty(res.bad_probe_idx)
        deconv_type(res.bad_probe_idx,:) = [];
    end    
    deconv_type = cat(1,deconv_type(:));
    %split into the data of the types that you want    
    if strcmp(statname,'rho')
%        temp{i} = fisherZ([res.deconv_stats(strcmp(deconv_type,comptype)).rho]); %rho
       temp{i} = [res.deconv_stats(strcmp(deconv_type,comptype)).rho]; %rho-z
    elseif strcmp(statname,'err')
       temp{i} = [res.deconv_stats(strcmp(deconv_type,comptype)).err]; %mse
    elseif strcmp(statname,'s')
       temp{i} = [res.deconv_stats(strcmp(deconv_type,comptype)).s]; %skew 
    else
       error('unknown statname');
    end   
end %loop
%equate their sizes
maxlen = max(cellfun(@(x) size(x,2),temp,'UniformOutput',1));
temp = cellfun(@(x) [x, NaN(size(x,1),maxlen-size(x,2))],temp,'UniformOutput',0);
data = cat(1,temp{:});
end %function

%% gather data but split by methods
function [data,lags,xcorr_trace,type] = LoadResults(fn,deconv_type)
%each cell of data is a different measure, column within cell is deconv
%type
res = load(fn);
if nargin <2
   deconv_type = {'lr_gcamp','glm','feedforward','none'}; 
end
deconv_type = repmat(deconv_type,size(res.test_idx,1),1);
if ~isempty(res.bad_probe_idx)
    deconv_type(res.bad_probe_idx,:) = [];
end
deconv_type = cat(1,deconv_type(:));
%split into the data of the types that you want
type = {'none','lr_gcamp','glm','feedforward'}; %types to use and order to plot
data = cell(1,3);
xcorr_trace = cell(1,numel(type));
lags = NaN(sum(strcmp(deconv_type,type{1})),numel(type));
for i = 1:numel(type)
%    data{1}(:,i) = fisherZ([res.deconv_stats(strcmp(deconv_type,type{i})).rho]); %rho-z
   data{1}(:,i) = [res.deconv_stats(strcmp(deconv_type,type{i})).rho]; %rho
   data{2}(:,i) = [res.deconv_stats(strcmp(deconv_type,type{i})).err]; %mse   
   data{3}(:,i) = ([res.deconv_stats(strcmp(deconv_type,type{i})).s]); %absolute skew
   lags(:,i) = [res.deconv_stats(strcmp(deconv_type,type{i})).x_lag]; %lag
   xcorr_trace{i} = [res.deconv_stats(strcmp(deconv_type,type{i})).xcorrvect]; %xcorr trace
end
xcorr_trace{i+1} = [res.deconv_stats(strcmp(deconv_type,type{i})).autocortrue];
end %function
%%
function [data,lags,xcorr_trace,type] = LoadResultsTrain(fn,deconv_type)
%each cell of data is a different measure, column within cell is deconv
%type
res = load(fn);
if nargin <2
   deconv_type = {'lr_gcamp','glm','feedforward','none'}; 
end
deconv_type = repmat(deconv_type,size(res.test_idx,1),1);
if ~isempty(res.bad_probe_idx)
    deconv_type(res.bad_probe_idx,:) = [];
end
deconv_type = cat(1,deconv_type(:));
%split into the data of the types that you want
type = {'none','lr_gcamp','glm','feedforward'}; %types to use and order to plot
data = cell(1,3);
xcorr_trace = cell(1,numel(type));
lags = NaN(sum(strcmp(deconv_type,type{1})),numel(type));
for i = 1:numel(type)
%    data{1}(:,i) = fisherZ([res.deconv_stats(strcmp(deconv_type,type{i})).rho]); %rho-z
   data{1}(:,i) = [res.deconv_stats_train(strcmp(deconv_type,type{i})).rho]; %rho
   data{2}(:,i) = [res.deconv_stats_train(strcmp(deconv_type,type{i})).err]; %mse   
   data{3}(:,i) = ([res.deconv_stats_train(strcmp(deconv_type,type{i})).s]); %absolute skew
   lags(:,i) = [res.deconv_stats_train(strcmp(deconv_type,type{i})).x_lag]; %lag
   xcorr_trace{i} = [res.deconv_stats_train(strcmp(deconv_type,type{i})).xcorrvect]; %xcorr trace
end
xcorr_trace{i+1} = [res.deconv_stats_train(strcmp(deconv_type,type{i})).autocortrue];
end %function

%% make base plots
function handles = makeViolinPlots(data,spreadflag,fp,sz,col,labels)    
    titlestr = {'Correlation','MSE','Skew'};
    ylabelstr = {'Rho-z','MSE','\Delta Skew'};
    for i = 1:numel(titlestr)
        figure; hold on;
        vp = CompareViolins(data{i}',fp,'col',col,'connectline',[0.25 0.25 0.25 0.50],'plotspread',spreadflag);
        ylabel(ylabelstr{i}); xlabel('Method')
        set(gca,'xticklabels',labels)
        fp.FormatAxes(gca); grid on
        fp.SetTitle(gca,titlestr{i})
        fp.FigureSizing(gcf,sz,[])
        if i==2 && size(data{i},2)>2 %for the full plots
            yval = get(gca,'ylim');
            ylim([0 yval(2)]);
        elseif i==3
            xval = get(gca,'xlim');
            plot(xval,[0 0],'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--')
        end
    end %figure loop
    handles = get(groot, 'Children');
end
%% xcorr plots
function handles = makeXcorrPlots(xcorr_trace,fp,sz,col,labels)
    st = xcorr_trace{end};
    xcorr_trace=xcorr_trace(1:end-1);
    %equal axes for large plot
    ymin = min(cellfun(@(x) min(nanmean(x,2)-sem(x,2)), xcorr_trace,'UniformOutput',1));
    ymax = max(cellfun(@(x) max(nanmean(x,2)+sem(x,2)), xcorr_trace,'UniformOutput',1));

    %Camden - can add the true autocorrelation (yyaxis) for reference
    for i = 1:numel(xcorr_trace)
        figure; hold on; 
        y = nanmean(xcorr_trace{i},2);
        yt = nanmean(st,2);
        x = -1*(numel(y)-1)/2:(numel(y)-1)/2;
        yerr = sem(xcorr_trace{i},2);
        %get center samples
        idx = ismember(x,-150:150);
        x = x(idx==1);
        yt = yt(idx==1);
        y = y(idx==1);
        yerr = yerr(idx==1);   
        plot(x,yt,'linewidth',fp.line_width,'color',[0.1 0.1 0.1 0.5]);
        shadedErrorBar(x,y,yerr,'lineprops',{'linestyle','-','color',col{i}});    
        ylim([ymin,ymax]);
        yval = get(gca,'ylim');
        plot([0 0],yval,'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--')
        plot([x(1) x(end)],[0 0],'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--')
        ylabel('Rho')
        xlabel('Samples');
        fp.FormatAxes(gca); grid on
        fp.SetTitle(gca,labels{i})
        fp.FigureSizing(gcf,sz{1},[])  
        xlim([-150,150])
%         yyaxis right
%         plot(x,yt,'linewidth',fp.line_width,'color',[0.1 0.1 0.1 0.5]);        

        %zoomed in
        figure; hold on; 
        y = nanmean(xcorr_trace{i},2);   
        yt = nanmean(st,2);
        x = -1*(numel(y)-1)/2:(numel(y)-1)/2;    
        yerr = sem(xcorr_trace{i},2);
        %get center 10 samples
        idx = ismember(x,-5:5);
        x = x(idx==1);
        yt = yt(idx==1);
        y = y(idx==1);
        yerr = yerr(idx==1);      
%         plot(x,yt,'linewidth',fp.line_width,'color',[0.1 0.1 0.1 0.5]);
        shadedErrorBar(x,y,yerr,'lineprops',{'linestyle','-','color',col{i}});        
        yval = get(gca,'ylim');
        plot([0 0],yval,'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--')
        plot([x(1) x(end)],[0 0],'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--')
        ylabel('Rho')
        xlabel('Samples');
        fp.FormatAxes(gca); grid on
        fp.SetTitle(gca,labels{i})
        fp.FigureSizing(gcf,sz{2},[])   
        ylim(yval);
%         yyaxis right
%         p=plot(x,yt,'linewidth',fp.line_width,'color',[0.1 0.1 0.1 0.5]);
%         set(gca,'YColor','k');
    end
    handles = get(groot, 'Children');
end
%% statistical comparisons
function handles = makeStatsPlots(data,fp,sz,labels)    
    titlestr = {'Correlation','MSE','Skew'};
    for i = 1:numel(titlestr)
        n=size(data{1},2);
        figure; hold on; 
        idx = nchoosek(1:n,2);
        p_mat = NaN(size(data{i},2));
        for j = 1:size(idx,1)   
           p_mat(idx(j,1),idx(j,2)) = -log10(signrank(data{i}(:,idx(j,1)),data{i}(:,idx(j,2))));
        end
        imagesc(p_mat,[0 4])

        colormap(flipud(gray(4))); 
        c=colorbar('Ticks',0:1:4);
        c.Label.String = '-log_1_0(p)';
        c.Label.Interpreter = 'tex';        
        box on; 
        set(gca,'Xtick',1:numel(labels),'Ytick',1:numel(labels),'XTickLabel',...
            labels,'YTickLabel',labels,'XTickLabelRotation',90)
        fp.FormatAxes(gca);
        xlim([0.5 n+0.5]); ylim([0.5 n+0.5])
        %add grid lines
        for j = 1:4
          line([0, n+0.5], [j+0.5, j+0.5], 'Color', 'w','linewidth',2);
          line([j+0.5, j+0.5],[0, n+0.5], 'Color', 'w','linewidth',2);
        end 
        fp.SetTitle(gca,titlestr{i})
        fp.FigureSizing(gcf,sz,[])
    end %figure loop
    handles = get(groot, 'Children');
end


