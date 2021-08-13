%% Deconvolution Figure Pipeline
%Camden MacDowell - timeless


%% Figure 2
fp = fig_params_deconvolutionpaper; 
col = [{fp.c_none},{fp.c_lr},{fp.c_glm},{fp.c_ff}];
%load example traces and plot

%load fit to 600uM 
fn = GrabFiles('within_compare\w*std\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution'});
[data,lags,xcorr_trace,type] = LoadResults(fn{end-1});%get the desired depth 
labels = {'None','LR','GLM','fNN'};

%make plots for rho, err, and skew
titlestr = {'Correlation','MSE','Skew'};
ylabelstr = {'Rho-z','MSE','\Delta Skew'};
for i = 1:numel(titlestr)
    figure; hold on;
    vp = CompareViolins(data{i}',fp,'col',col,'connectline',[0.25 0.25 0.25 0.50],'plotspread',1);
    ylabel(ylabelstr{i}); xlabel('Method')
    set(gca,'xticklabels',labels)
    fp.FormatAxes(gca); grid on
    fp.SetTitle(gca,titlestr{i})
    fp.FigureSizing(gcf,[2 2 6 8],[])
    if i==2
        yval = get(gca,'ylim');
        ylim([0 yval(2)]);
    end
end %figure loop

%plot the 
%also could look at the width of the xcorrelation as an idea of how well it
%matches the spiking dynamics. 



%%
function [data,lags,xcorr_trace,type] = LoadResults(fn)
%each cell of data is a different measure, column within cell is deconv
%type
res = load(fn);
deconv_type = {'narx','lr_gcamp','lr_glm','glm','feedforward','none'}; 
deconv_type = repmat(deconv_type,size(res.test_idx,1),1);
deconv_type = cat(1,deconv_type(:));
%split into the data of the types that you want
type = {'none','lr_gcamp','glm','feedforward'}; %types to use and order to plot
data = cell(1,3);
xcorr_trace = cell(1,numel(type));
lags = NaN(sum(strcmp(deconv_type,type{1})),numel(type));
for i = 1:numel(type)
   data{1}(:,i) = fisherZ([res.deconv_stats(strcmp(deconv_type,type{i})).rho]); %rho
   data{2}(:,i) = [res.deconv_stats(strcmp(deconv_type,type{i})).err]; %mse   
   data{3}(:,i) = ([res.deconv_stats(strcmp(deconv_type,type{i})).s]); %absolute skew
   lags(:,i) = [res.deconv_stats(strcmp(deconv_type,type{i})).x_lag]; %lag
   xcorr_trace{i} = [res.deconv_stats(strcmp(deconv_type,type{i})).xcorrvect]; %xcorr trace
end
end %function