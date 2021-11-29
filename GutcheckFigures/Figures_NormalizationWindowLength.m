function Figures_NormalizationWindowLength() 
%Camden MacDowell - timeless
fp = fig_params_deconvolutionpaper;
%load the data
fn = GrabFiles('\w*block\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Normalization\WindowSweepResults'});

%get the duration
dur = cellfun(@(x) load(x,'dur'),fn,'UniformOutput',1); 
dur = [dur(:).dur];
unique_dur = unique(dur);
%loop through each method and plot the training results
train = cellfun(@(x) load(x,'deconv_stats_train'),fn,'UniformOutput',0); 
train = cellfun(@(x) [x.deconv_stats_train(:).rho],train,'UniformOutput',0);
train = cat(1,train{:});
test = cellfun(@(x) load(x,'deconv_stats'),fn,'UniformOutput',0); 
test = cellfun(@(x) [x.deconv_stats(:).rho],test,'UniformOutput',0);
test = cat(1,test{:});

type = load(fn{1},'type'); type = type.type;
n_probe = 4; 
type_idx = repmat(1:numel(type),n_probe,1);
type_idx = type_idx(:);

data_dur = {};
for i = 1:numel(type)
   %organize by method and duration 
   temp_mat = {};
   for j = 1:numel(unique_dur)
      temp = train(dur==unique_dur(j),type_idx==i);
      temp_mat{j} = temp(:);
   end   
   data_dur{i} = [temp_mat{:}];
end

figure; hold on; 
col = {fp.c_lr,fp.c_glm,fp.c_ff,fp.c_none};
for i = 1:numel(data_dur) %loop through each method
   n = size(data_dur{i},2);
   shadedErrorBar(1:n,nanmean(data_dur{i}),sem(data_dur{i}),'lineprops',{'color',col{i},'linewidth',2});      
end
set(gca,'xtick',1:n,'XTickLabel',unique_dur,'xlim',[0.5,numel(unique_dur)+0.5]);
ylabel('Rho')
xlabel('Window Duration (sec)')
fp.FigureSizing(gcf,[3 2 6 6],[])
fp.FormatAxes(gca);
p.LineWidth=1; box on 
saveCurFigs(gcf,{'-dpng','-dsvg'},'normalizationwindow','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Normalization',0); close all

%% out of curiousity, plot the glm kernels
kernel = cellfun(@(x) load(x,'trained_opts'),fn,'UniformOutput',0); 
kernel = cellfun(@(x) [x.trained_opts{:}.glmkernel],kernel,'UniformOutput',0);

%%
%plot the kernels for each distance
col = repmat(linspace(0,0.7,numel(unique_dur)),3,1);
figure; hold on; 
for i = 1:numel(unique_dur)    
    temp = [kernel{dur==unique_dur(i)}];
    shadedErrorBar(1:size(temp,1),nanmean(temp,2),sem(temp,2),'lineprops',{'color',col(:,i),'linewidth',2});                
end %function end
saveCurFigs(gcf,{'-dpng','-dsvg'},'normalizationwindow_glmkernel','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Normalization',0); close all
%%
    
    



























