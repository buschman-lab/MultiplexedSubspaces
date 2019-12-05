%Saving Information
savedir = 'C:\Users\macdo\OneDrive\Buschman Lab\AnalysisCode_Repository\Mesoscale Network Dynamics 2019 Analyses\Training and Testing Statistics\11-5-2019';
if ~exist(savedir)
    mkdir(savedir)
end
savefigs = 1;
writestats = 1;  
filename = [savedir 'FigureStatistics.csv'];
%Delete the data table and rewrite
if writestats || exist(filename)
    delete(filename);
end    
base = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\';
fp = fig_params; 

mouse_num = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO\MouseNumInfoByRec.mat');
group = isVPA(mouse_num.mousenum); 

%% Figure comparing cnmf to nmf and pca
other_methods = CompileStats(GrabFiles('block',0,{[base 'TrainRepitoires\TrainingFit_CompareDiscoveryMethods_full']}),{'nmf','spca'},0,group);
cnmf = CompileStats(GrabFiles('block',0,{[base 'TrainRepitoires\TrainingFit_Lambda4e-4_Kval28']}),{'ExpVar_all','numFactors'},0,group);
%%
stats = Plot_CompareDiscoveryMethods(other_methods,cnmf); 

if writestats 
    save([savedir 'comparediscovermethodsstats.mat'],'stats');
end
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','Train_ExpVar',savedir,1);
   close all
end

%% Distribution of Events
D = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28';
params = gatherHs(318,'TrainRepitoire_block',0,D);
%Get the average number for all factors in a epoch
%%
figure('position',[680   451   268   527]); hold on;
setFigureDefaults();
occurance = [params(:).IndiOccurance]/2;
hh = histogram(occurance,ceil(max(occurance)),'BinEdges',(0:1:ceil(max(occurance))));
set(hh(1),'FaceColor',fp.c_discovery,'FaceAlpha',0.4,'EdgeColor',fp.c_discovery,'EdgeAlpha',0.8);
ymax = get(gca,'ylim');
line([mean(occurance),mean(occurance)],[0,ymax(2)+2],'LineWidth',2,...
    'linestyle','--','Color',[0.3 0.3 0.3]);
ylabel('# Motifs');
xlabel({'Motif Frequency';'(occurence/minute)'})
ylim([0 ymax(2)+2]);
xlim([0 max(occurance)])
set(gca,'position',[2 3 3.5 8.5])
title({'Motif Frequency';''},'FontWeight','normal','units','centimeters','position',[1.75,9.25],'FontName','Arial');
text(mean([params(:).Occurance])+1,400,['\mu=', num2str(round(mean(occurance),2))],'FontSize',16,'FontWeight','normal','FontName','Arial');

%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','FittingStats_NumEvents',savedir,1);
   close all
end
if writestats
   u = mean(occurance);
   ci = bootci(1000,@mean,occurance);
   temp = sort(occurance,'ascend');   
   temp = temp(floor(numel(temp)/2)-1);
   T = {'Variable','MotifOccurance_Fig1C';'Mean',u;'CI_l',ci(1);'CI_u',ci(2);'HalfOccurAtLeastPerMin',temp};
   [T_new] = WriteStatToTable(T,filename,0);       
end

%% Figure 2 D Loadings
%organize data
data = CompileStats(GrabFiles('block',0,{[base 'TrainRepitoires\TrainingFit_Lambda4e-4_Kval28']}),{'ExpVarLoad'},0,group);
data = {data(:).ExpVarLoad};
data = MakeCellsEqual(data);
data = cumsum(data*100,1,'omitnan');
ci = NaN(size(data,1),2);
for i = 1:size(data,1)
    ci(i,:) = bootci(1000,@nanmedian,data(i,:));
end

%%
figure('position',[1092 252 336 586]); hold on
errorbar(1:1:size(data,1),nanmedian(data,2),(nanmedian(data,2)-ci(:,1)),(ci(:,2)-nanmedian(data,2)),'LineWidth',2,'Color',[0.4 0.4 0.4])
title({'All Motifs Are Used to';'Explain Neural Activity'},'FontWeight','normal','units','centimeters','position',[2.5,9.25]);
xlabel({'Motif';'(Ordered By Decreasing';'Variance Per Epoch)'})
ylabel({'Relative Percent Explained';'Variance (Cummulative)'});
ylim([0 100]);
xticks((0:5:size(data,1)))
xlim([0 size(data,1)+0.5])
setFigureDefaults();
set(gca,'ytick',(0:20:100));
set(gca,'position',[3 4 5 8.5]);
%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','FittingStats_FactorLoadings',savedir,1);
   close all
end
if writestats
   for i = 1:size(data,1)
       u = nanmedian(data(i,:));       
       n = sum(~isnan(data(i,:)));
       T = {'Variable',sprintf('Loading%d_Fig1D',i);'Mean',u;'CI_l',ci(i,1);'CI_u',ci(i,2);'NumDataPoints',n};
       [T_new] = WriteStatToTable(T,filename,0);
   end
end
%% Motifs explained variance within and between animals
dir_str = {'TestRepitoires\SameAnimal_Fit_Kval28_lambda1\','TestRepitoires\DifferentAnimal_Fit_Kval28_Lambda1'};
data = [];
for i = 1:numel(dir_str)
    temp = CompileStats(GrabFiles('_ModFlag0_',0,{[base dir_str{i}]}),{'ExpVar_all'},0,group);
    if i == 2%get average of the multiple runs for different animals
        temp = nanmedian(cell2mat(cat(1,temp(:).ExpVar_all)),2);
    else        
        temp = cell2mat([temp(:).ExpVar_all]);
    end
    data(i,:) = temp*100; 
end   
%%
figure('position',[680  200  242   600]); hold on
vp =CompareViolins(data,fp,'label',{'Same Animal','Different Animal'},'col',{[0 0 0.75],[0.3 0.7 .4]});
set(gca,'XTickLabelRotation',45)
ylim([0 100]);
ylabel({'Percent Explained Variance'})
setFigureDefaults;
set(gca,'position',[2,4,4,8.5])
[p, h] = signrank(data(1,:),data(2,:));
AddSig(h,p,[1,2,100,100],2,5,1)

if writestats
    clear stats;
    stats.within_median = median(data(1,:));
    stats.ci_within = bootci(1000,@median,(data(1,:)));
    stats.between_median = median(data(2,:));
    stats.ci_between = bootci(1000,@median,(data(2,:)));
    stats.difference  =median(data(1,:))-median(data(2,:));
    stats.pval = signrank(data(1,:),data(2,:));
    save([savedir filesep 'within_between_stats.mat'],'stats');
end
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','Within_between_expvar',savedir,1);
   close all
end


%% End Figure 2





%% Figure 3: Dynamics
%Load All motifs
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\ClusteredDPs_rawdata.mat','W_all');
label={'Static','SpatialStatic','framewise'};
%Smooth everything
[nP,nK,nT] = size(W_all);
W_all_smooth = W_all;
for cur_k = 1:nK
    temp = squeeze(W_all(:,cur_k,:));
    temp = reshape(temp,sqrt(nP),sqrt(nP),nT);
    temp = SpatialGaussian(temp);
    W_all_smooth(:,cur_k,:) = reshape(temp,nP,1,nT);
end
Comparison = toggleModFlag(W_all_smooth,[],1); %Completely static motifs
data = NaN(size(Comparison,2),2);
[data(:,1), data(:,2),~] = CorrelateMotifs(W_all_smooth,Comparison);
data = 1-data; %Dissimilarity 

figure('position',[680   514   242   464]); hold on
COL = {[0.4 0.4 0.4]};
vp =CompareViolins(data(:,2)',fp,'col',COL,'label','');
ylabel({'Dissimilarity (1-\rho)';'With Dynamic Motifs'})
xlabel({'Static';'Network'});
setFigureDefaults();
set(gca,'position',[3 3 2 8.5],'ylim',[0 1])

clear W_all W_all_smooth
%%

if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','DynamicVsStaticComparison_Figure2',savedir,1);
   close all
end
if writestats
   for i = 1:size(data,2)
       u = median(data(:,i));
       ci = bootci(1000,@median,data(:,i));
       T = {'Variable',sprintf('Dissimilarity_%d_Fig1E',i);'Mean',u;'CI_l',ci(1);'CI_u',ci(2)};
       [T_new] = WriteStatToTable(T,filename,0);
   end
end  

%% Figure 2G temporal-Autocorrelation
cd('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28');
load('ClusteredDPs_rawdata.mat','W_all')
[metric,tau_hl,rmv_metric,rmv_tau] = DynamicMetric(W_all);

%compute the half life of the mean 
x = (1:1:(size(metric,2))*1)';
y = nanmedian(metric,1)';
f = fit(x,y,'exp1','StartPoint',[0,0]);
tau_hl_mean = -1*(log(2)/f.b);
    
tau_hl_bootstrap = zeros(1000,1);
for i = 1:1000
    temp = nanmedian(metric(randi(size(metric,1),size(metric,1),1),:),1);
    x_temp = (1:75:(numel(temp))*75)';
    y_temp = temp';
    f_temp = fit(x_temp,y_temp,'exp1','StartPoint',[0,0]);
    tau_hl_bootstrap(i) = -1*(log(2)/f_temp.b);
end
%%
figure('position',[542   420   336   555]); hold on;
data = metric(:,1:end)'; 
shadedErrorBar(1:1:size(data,1),nanmedian(data,2),(nanstd(data,[],2)),'lineprops',{'-','color',...
[0.4 0.4 0.4],'linewidth',2},'transparent',1,'patchSaturation',0.2);
p1 = plot(f,x,y,'k');
set(p1,'linewidth',2);
legend off;
title({'Correlation of Spatial';'Patterns in Motifs';'Decays Over Time'},...
    'FontWeight','normal','units','centimeters','position',[2.25,9]);
xlabel({'Temporal Offset';'Within Motif (ms)'})
ylabel('Correlation');
ylim([-0.2 0.62]);
xlim([0.5 size(data,1)+0.5])
ylim([-0.15 0.6])
setFigureDefaults();
set(gca,'XTickLabel',(1:2:13)*75,'XTick',1:2:13,'XTickLabelRotation',45)
set(gca,'Position',[2.25 3.05 5 8.5],'YColor','k')

pval_store = ones(size(data,1),1);
h= ones(size(data,1),1);
for i = 1:size(data,1)      
    [pval_store(i),h(i)] = signrank(data(i,:),0,'tail','right','alpha',0.01/size(data,1));
end
%Add pvalue line
line([find(h==1,1,'first'),find(h==1,1,'last')],[0.6 0.6],'linewidth',2,'color',[0.5 0.5 0.5]);
% text(0.1,-0.05,sprintf('T = %.2f +/- %.2fms',tau_hl_mean*75,std(tau_hl_bootstrap)*75),'fontsize',16,'FontName','Arial'); %huge tail so use median

[pval, ~] = signrank((data(1,:)),(data(end,:)));
%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','Dynamicness',savedir,1);
   close all
end

if writestats
   fileID = fopen([savedir 'FigureStatistics_AutoCorrelation.txt'],'w');
   fprintf(fileID,'\nTemporalAutoCorrelation\n') 
   fprintf(fileID,'\nmeans = %.4g',...
       nanmedian(data,2));
   fprintf(fileID,'\nTemporalAutoCorrelation\n STD = %.4g',...
       nanstd(data,[],2));
   fprintf(fileID,'\nTemporalAutoCorrelation\n pval from 0')
   fprintf(fileID,'\npval = %.4g',...
       pval_store);   
   temp = data(1,:) - data(end,:);   
   ci = bootci(1000,@nanmedian,temp);    
   fprintf(fileID,'\ndifference between 2 and 12 = %.4g ci %d, %d, pval = %.4g',nanmedian(temp),ci(1),ci(2),pval);       
   fprintf(fildID,'\nFirst correlation %g ci %g and %g pval %d',nanmedian(data(1,:)),bootci(1000,@nanmedian,data(1,:)),signrank(data(1,:),0,'tail','right'))
   fprintf(fileID,'\nFirst correlation %g ci %g and %g pval %d',nanmedian(data(12,:)),bootci(1000,@nanmedian,data(12,:)),signrank(data(12,:),0,'tail','right'))   
   fprintf(fileID,'\nremoved tau (<4) = %g removed metric (no active) = %d',rmv_tau,rmv_metric);
   fprintf(fileID,'\nHalf-life of decay =%g 25th = %g, 75th',nanmedian(tau_hl),prctile(tau_hl,25),prctile(tau_hl,75))
   fprintf(fileID,'\nHalf-life of average decay =%g +/- SEM %d',tau_hl_mean*75,std(tau_hl_bootstrap),prctile(tau_hl,75))
   fclose(fileID); 
end
%% Static network vs motifs
rec_str = {'_ModFlag0_','_ModFlag1_'};
data = [];
for i = 1:numel(rec_str)
    temp = CompileStats(GrabFiles(rec_str{i},0,{[base 'TestRepitoires\SameAnimal_Fit_Kval28_lambda1\']}),{'ExpVar_all'},0,group);    
    temp = cell2mat([temp(:).ExpVar_all]);
    data(i,:) = temp*100; 
end   
%%
figure('position',[680  200  242   600]); hold on
vp =CompareViolins(data,fp,'label',{'Motif','Static Network'},'col',{[0 0 0.75],[0.5 0.5 0.5]});
set(gca,'XTickLabelRotation',45)
ylim([0 100]);
ylabel({'Percent Explained Variance'})
setFigureDefaults;
set(gca,'position',[2,4,4,8.5])
[p, h] = signrank(data(1,:),data(2,:));
AddSig(h,p,[1,2,100,100],2,5,1)

%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','Motif_vs_static_Network',savedir,1);
   close all
end


if writestats
    clear stats;
    stats.within_median = median(data(1,:));
    stats.ci_within = bootci(1000,@median,(data(1,:)));
    stats.between_median = median(data(2,:));
    stats.ci_between = bootci(1000,@median,(data(2,:)));
    stats.difference  =median(data(1,:))-median(data(2,:));
    stats.pval = signrank(data(1,:),data(2,:));
    save([savedir filesep 'within_static_dynamic.mat'],'stats');
end

%% END OF FIGURE 3 DYNAMICS










%% Basis Motifs
D = [];
D{1} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\SameAnimal_Fit_Kval28_Lambda1';
D{2} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\DifferentAnimal_Fit_Kval28_Lambda1';
D{3} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Kval28_Lambda1';
D{4} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Kval28_Lambda1';
base = {'Fitted_block_ModFlag0','Fitted_block_ModFlag0','Fitted_block_ModFlag0','Fitted_block_ModFlag1'};
label = {'Within Mice','Between Mice','Basis Motifs','Static Networks'};
data = [];
for i =1:numel(D)
    temp = gatherData(318,base{i},0,D{i},{'ExpVar_all'});
    temp = cell2mat(cat(1,temp(:).ExpVar_all));
    if size(temp,2)>1
       temp = nanmedian(temp,2);
    end
    data(i,:) = temp*100;
end
%%
figure('position',[150  150   550   550]); hold on;
COL = {[0 0 0.75],[0.25 0.75 0],[0.85,0.37,0.001],[0.5 0.5 0.5]};
vp =CompareViolins(data,fp,'label',label,'col',COL);
t=title({'Basis Motifs Capture';'Majority of Cortical Activity'},'Units','Centimeters',...
    'FontName','Arial','FontWeight','normal','VerticalAlignment','bottom','position',[4.2598 9 0]);
ylabel({'Percent Explained Variance'});
ylim([0 100])
setFigureDefaults();
set(gca,'XTickLabelRotation',45,'ytick',(0:20:100),'Clipping','off','position',[3 4 8.5 8.5])
[pval, h] = signrank(data(1,:),data(3,:));
AddSig(h,pval,[1 3 96.5],3,5,1); 
[pval, h] = signrank(data(1,:),data(2,:));
AddSig(h,pval,[1 2 24],3,6,1); 
[pval, h] = signrank(data(2,:),data(3,:));
AddSig(h,pval,[2 3 38],3,6,1); 
[pval, h] = signrank(data(3,:),data(4,:));
AddSig(h,pval,[3 4 12],3,6,1); 
%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-png','ClusteredFittingExplainedVariance',savedir,1);
   close all
end
if writestats    
    clear stats;
    stats.basis_median = median(data(3,:));
    stats.ci_basis = bootci(1000,@median,(data(3,:)));
    stats.staticbasis_median = median(data(4,:));
    stats.ci_staticbasis = bootci(1000,@median,(data(4,:)));
    stats.difference_staticbasistobasis  =median(data(3,:))-median(data(4,:));
    stats.pval_staticbasistobasis = signrank(data(3,:),data(4,:));
    
    stats.difference_orig_between =median(data(1,:))-median(data(2,:));
    stats.pval_difference_orig_between= signrank(data(1,:),data(2,:));
    
    stats.difference_orig_basis  =median(data(1,:))-median(data(3,:));
    stats.pval_orig_basis  = signrank(data(1,:),data(3,:));
    
    stats.difference_between_basis =median(data(2,:))-median(data(3,:));
    stats.pval_between_basis = signrank(data(2,:),data(3,:));
    save([savedir filesep 'basismotifs_stats.mat'],'stats');

end

%% Figure 3D Exp Var Loadings 
D = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Kval28_lambda1';
base = 'Fitted_block_ModFlag0';
params = gatherData(318,base,0,D,{'ExpVarLoad'});
data = {params(:).ExpVarLoad};
maxlength = max(cellfun(@numel, data));
data= cellfun(@(v) [v, nan(1, maxlength-numel(v))], data, 'UniformOutput', false);
data = vertcat(data{:})';
data(data==0)=NaN;
[~, idx] = sort(nanmedian(data,2),'descend'); %reorder by explained variance
data = data(idx,:)*100;
data = cumsum(data,1,'omitnan');
ci = NaN(size(data,1),2);
for i = 1:size(data,1)
    ci(i,:) = bootci(1000,@nanmedian,data(i,:));
end
%%
figure('position',[680   558   422   420]); hold on;
errorbar(1:1:size(data,1),nanmean(data,2),(nanmean(data,2)-ci(:,1)),(ci(:,2)-nanmean(data,2)),'LineWidth',2,'Color',[0.4 0.4 0.4])
title({'All Basis Motifs Are';'Used to Explain Activity'},'FontWeight','normal','FontName','Arial',...
    'units','centimeters','position',[2 5.5]); 
xlabel('Basis Motif');
ylabel({'Relative Percent';'Explained Variance'});
setFigureDefaults();
ylim([0 100]);
xticks((1:4:size(data,1)))
xlim([0 size(data,1)+0.5])
box off
set(gca,'Position',[4 2.25 4.5 5],'Clipping','off');

%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','BasisMotifLoadings',savedir,1);
   close all
end
if writestats
   for i = 1:size(data,1)
       u = nanmean(data(i,:));       
       n = sum(~isnan(data(i,:)));
       T = {'Variable',sprintf('Loadings%d_Fig3',i);'Mean',u;'CI_l',ci(i,1);'CI_u',ci(i,2);'NumDataPoints',n};
       [T_new] = WriteStatToTable(T,filename,0);
   end
end

%% Kvalues for Basis Motifs
D = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\XcorrClusterings';
COL = [0.2 0.2 0.2];
[knnval, kval, numClust] = AnalyzeKSweep(D,15);
klist = unique(knnval);
%%
figure(); hold on
for i = 1:length(klist)
    plot(kval(knnval==klist(i)),numClust(knnval==klist(i)),'Marker','.','color',COL(i,:),'MarkerSize',20,'LineStyle','none')
end
title({'Basis Motifs Are Robust';'to CNMF Parameters'},'FontWeight','normal','FontName','Arial',...
    'units','Centimeters','Position',[2 5.5])
xlim([0 40])
ylim([1 20])
xlabel({'Max # Discoverable';'Motifs Per Epoch ({\itK})'});
ylabel('# Basis Motifs');
set(gca,'LineWidth',2,'FontSize',16,'FontWeight','normal','FontName','Arial')
setFigureDefaults;
set(gca,'units','centimeters','position',[4 3 4.5 5]);
plot(28,numClust(kval==28),'Marker','o','color',[0.9 0 0],'Markersize',15,'LineStyle','none','Linewidth',2);
set(gcf,'Position',[680   438  500  375]);
%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','KSweep',savedir,1);
   close all
end
%% END BASIS MOTIF FIGURE










%% SOCIAL FIGURE
D = [];
D{1} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Kval28_lambda1';
D{2} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\Social\Social_TestRepitoires\NewClustered_Fit_Kval28_lambda1';
D{3} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\Social\Social_TestRepitoires\NewsClustToRegular_Fit_Kval28_lambda1';
D{4} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\Social\Social_TestRepitoires\NewsClustToSocial_Fit_Kval28_lambda1/';

base = {'Fitted_block_ModFlag0','Fitted_block_ModFlag0','Fitted_block_ModFlag0','Fitted_block_ModFlag0'};
label = {'Solo: Solo','Solo: Social','Social: Solo','Social: Social'};

params = cell(1,numel(D));

for i =1:numel(D)
    if i == 1
        params{i} = gatherData(318,base{i},0,D{i},{'ExpVar_all'}); %original --> original
    elseif i == 2
        params{i} = gatherData(123,base{i},2,D{i},{'ExpVar_all'}); %original --> social
    elseif i == 3
        params{i} = gatherData(318,base{i},0,D{i},{'ExpVar_all'}); %social --> Original
    elseif i == 4
        params{i} = gatherData(123,base{i},2,D{i},{'ExpVar_all'}); %social --> Social
    end
end
%%
% Compare Explained Variance
figure; hold on
data = NaN(numel(D),144); %Have to nanpad since unequal sizes
for i = 1:numel(D)
    temp = [params{i}];
    temp = cell2mat([temp(:).ExpVar_all]);
    data(i,1:numel(temp)) = temp*100; 
end
COL = {[0.85,0.37,0.001],[0.85,0.37,0.001],[.14 0.4 0.6431],[0.14 0.4 0.6431]};
vp =CompareViolins(data,fp,'label',label,'col',COL);
title({'Basis Motifs Generalize Across';'Social and Solo Environments'},'FontWeight','normal','FontName','Arial',...
    'units','centimeters','Position',[4 9.25]);
ylabel({'Percent Explained Variance'});
ylim([0 100]);
setFigureDefaults;
set(gca,'XTickLabelRotation',45,'ytick',(0:20:100))
[pval, h] = ranksum(data(1,:),data(2,:)); %independent samples solo vs social environment
AddSig(h,pval,[1 2 99],2,5,1); 
[pval, h] = ranksum(data(3,:),data(4,:)); %Also independent
AddSig(h,pval,[3 4 99],2,5,1); 

xlabel({'Basis Motif Type: Environment'})
set(gca,'Clipping','off','LineWidth',2,'Fontsize',16,'Fontweight','normal');
%Make it so all Y axis are identical across figures
set(gca,'Units','centimeters','position',[3 5 8 8.5]);
set(gcf,'Position',[524   290   500  600]);
%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-png','SocialExpVar_BothDirections',savedir,1);
   close all
end
if writestats
   for i = 1:size(data,1)
       u = nanmedian(data(i,:));       
       ci = bootci(1000,@nanmedian,data(i,:));
       T = {'Variable',sprintf('SocialExpVar%d_Fig4',i);'Mean',u;'CI_l',ci(1);'CI_u',ci(2)};
       [T_new] = WriteStatToTable(T,filename,0);
   end   
   fileID = fopen([savedir 'FigureStatistics.txt'],'a');
   fprintf(fileID,'\nSocialMotifFittingPEV\n');
   temp = (data(1,:))-(data(2,:));
   ci = bootci(1000,@nanmedian,temp);
   fprintf(fileID,'\nSolo: Diff = %.4g CI %.4g %.4g, pval=%.4g',...
       nanmedian(temp),ci(1),ci(2), signrank(data(1,:),data(2,:)));
   temp = (data(3,:))-(data(4,:));
   ci = bootci(1000,@nanmedian,temp);
   fprintf(fileID,'\nSocial: Diff = %.4g CI %.4g %.4g, pval=%.4g',...
       nanmedian(temp),ci(1),ci(2), signrank(data(1,:),data(2,:)));
   fclose(fileID); 
end

%% Social relative explained variance vs Solo
D = [];
D{1} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Kval28_lambda1';
D{2} = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\Social\Social_TestRepitoires\NewClustered_Fit_Kval28_lambda1';
base = {'Fitted_block_ModFlag0','Fitted_block_ModFlag0'};
label = {'Individual','Social'};
params = cell(1,numel(D));
ExpVar = cell(1,numel(D));
for i =1:numel(D)
    if i == 2
        params = gatherData(123,base{i},2,D{i},{'ExpVarLoad'}); %Social
        data = {params(:).ExpVarLoad};
        maxlength = max(cellfun(@numel, data));
        data= cellfun(@(v) [v, nan(1, maxlength-numel(v))], data, 'UniformOutput', false);
        data = vertcat(data{:})'; 
        data = data(idx,:); %reorder according to exp var in original data set
        ExpVar{i} = data;
    else
        params = gatherData(318,base{i},0,D{i},{'ExpVarLoad'}); %original
        data = {params(:).ExpVarLoad};
        maxlength = max(cellfun(@numel, data));
        data= cellfun(@(v) [v, nan(1, maxlength-numel(v))], data, 'UniformOutput', false);
        data = vertcat(data{:})';
        data(data==0)=NaN;
        [~, idx] = sort(nanmedian(data,2),'descend'); %sort by the exp var in the original dataset
        data = data(idx,:); %reorder according to descendign occurance        
        ExpVar{i} = data;
    end
end
%%
figure(); hold on
%Plot the mean of both (pad with N to
plot(0:.1:1,0:.1:1,'LineWidth',1,'Color',[0.5 0.5 0.5],'LineStyle','--')
scatter(nanmedian(ExpVar{1},2),nanmedian(ExpVar{2},2),20,'filled','k')

%Add circles around the significant ones, were significance color
h_store = zeros(1,size(ExpVar{1},1));
for i = 1:size(ExpVar{1},1)
    [pval, h] = ranksum(ExpVar{1}(i,:),ExpVar{2}(i,:));
    if h == 1                
        if pval<0.0001
            h_store(i) = 1; %save off Hs for binomial prob calc
            scatter(nanmedian(ExpVar{1}(i,:)),nanmedian(ExpVar{2}(i,:)),100,[1 0 0],'filled')
        elseif pval<0.05/numel(idx)
            h_store(i) = 1; %save off Hs for binomial prob calc
            scatter(nanmedian(ExpVar{1}(i,:)),nanmedian(ExpVar{2}(i,:)),50,[1 0 0],'filled')          
        end
    else
        h_store(i) = 0; %save off Hs for binomial prob calc
    end
    text([nanmedian(ExpVar{1}(i,:))-.005,nanmedian(ExpVar{1}(i,:))-.005],...
    [nanmedian(ExpVar{2}(i,:)),nanmedian(ExpVar{2}(i,:))],...
    sprintf('%d',i),'Fontsize',16,'FontWeight','normal',...
    'Color',[0.3 0.3 0.3],'HorizontalAlignment','center') 
    
end

title({'Basis Motifs Are Expressed Differently';'in Solo and Social Environments'},'FontWeight','normal');
xlabel({'Relative Percent Explained';'Variance Solo Environment'});
ylabel({'Relative Percent Explained';'Variance Social Environment'});
ylim([0 0.26]);
xlim([0 0.26]);
setFigureDefaults;
set(gca,'xtick',get(gca,'xtick'),'xticklabels',get(gca,'XTick')*100,'ytick',get(gca,'YTick'),'yticklabels', get(gca,'YTick')*100);
set(gca,'Position',[3 3 8.5 8.5]);
set(gcf,'Position',[652   272   502   518])

%Binomial Test
pval_bin = 1 - binocdf(sum(h_store),numel(h_store),0.05/size(data,1));
%%
if writestats
   fileID = fopen([savedir 'FigureStatistics.txt'],'a');
   fprintf(fileID,'\nSocialVSSoloExpression\n num diff = %g, pval = %.4g',sum(h_store),pval_bin);
   fclose(fileID); 
end

if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','SimpleVsComplexExpression',savedir,1);
   close all
end
%% Figure 4 D Cummulative Explained Variance
% Plot the explained variance in descending order for both solo (expvar1) and social (expvar2)
figure('position',[ 706   192   520   562]); hold on; 
data = ExpVar{2};
[~, idx] = sort(nanmedian(data,2),'descend');
data = data(idx,:);

%Convert to cumsum
data1 = cumsum(ExpVar{1},1,'omitnan');
data2 = cumsum(data,1,'omitnan');

ci_data1 = [];
ci_data2 = [];
for i = 1:numel(idx)
    ci_data1(:,i) = bootci(1000,@nanmedian,data1(i,:));
    ci_data2(:,i) = bootci(1000,@nanmedian,data2(i,:));
end

%Simple
errorbar(1:1:numel(idx),nanmedian(data1,2),nanmedian(data1,2)-ci_data1(1,:)',ci_data1(2,:)'-nanmedian(data1,2),'LineWidth',1.5,'Color',[0.4 0.4 0.4])
%Social
errorbar(1:1:numel(idx),nanmedian(data2,2),nanmedian(data2,2)-ci_data2(1,:)',ci_data2(2,:)'-nanmedian(data2,2),'LineWidth',1.5,'Color',[1 0.4 0.6])

%Mark sig differences
solo_cs = cumsum(ExpVar{1},1,'omitnan');
social_cs = cumsum(data,1,'omitnan');
for i = 1:numel(idx)
    [pval h] = ranksum(solo_cs(i,:),social_cs(i,:));
    if h == 1 
        pos = min(nanmedian(solo_cs(i,:)),nanmedian(social_cs(i,:)));
        if pval<0.0001                        
            text([i i],[pos-.1,pos-.1],'**','Color',[0.1 0.1 0.1],'FontSize',16,'FontWeight','normal','HorizontalAlignment','Center');
        elseif pval<0.05/numel(idx)
            text([i i],[pos-.1,pos-.1],'*','Color',[0.1 0.1 0.1],'FontSize',16,'FontWeight','normal','HorizontalAlignment','Center');
        end
    end
end

legend('Solo','Social','Location','best')
title({'Fewer Basis Motifs Are';'Needed to Explain Neural';'Activity in Social Environment'},'FontWeight','normal');
xlabel({'Basis Motif';'(Ordered by Decreasing';'Percent Explained Variance)'});
ylabel({'Relative Percent Explained';'Variance (Cummulative)'})
ylim([0 1])
setFigureDefaults;
set(gca,'ytick',get(gca,'YTick'),'yticklabels', get(gca,'YTick')*100,'Xtick',[1,size(data,1)],'xticklabels',{'first','last'});
xlim([0 size(data,1)+0.5])
set(gca,'Units','centimeters','position',[3.6 3.5 8.5 8.5]);
%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','Suppliment_SimpleVsComplexOccuranceExpVar',savedir,1);
   close all
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    
function params = gatherHs(numIter,base,grpflag,directory)
    fprintf('\n Compiling data...\n');
    %Go to the directory and load the recordings, divide up by group and    
    cd(directory)
    load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO\MouseNumInfoByRec.mat')
    %Get grp info
    grp = isVPA(mousenum);
    if grpflag == 1 %do vpa animals
        rmvFlag = 0; %remove Sal
    elseif grpflag ==0 %do saline animals
        rmvFlag = 1; %remove VPA
    else
        rmvFlag =2;
    end
    params = [];
    for i = 1:numIter
        tempval = load(sprintf('%s_%d.mat',base,i),'H'); 
        params(i).H = tempval.H;
        occur = [];
        for cur_H = 1:size(tempval.H,1)
            temp = tempval.H(cur_H,:);
            temp(temp>(nanmedian(temp,2)+(1*std(temp,[],2)))) = 1;
            temp(temp~=1) = 0; 
            temp = diff(temp); 
            occur(cur_H) = numel(find(temp==1));
        end %h loop
        params(i).IndiOccurance = occur; 
        params(i).Occurance = mean(occur);
    end %iter loop
    

    params(grp==rmvFlag)=[];
end %function end

    
    
    
    
    
    
    
    
    