function Plot_MotifK()
% Camden MacDowell - timeless
% plot the one second motif duration, all the motif durations, and compares
% data is a 1 x N cell array, where each c

% optionally load data. Hardcoded for 2020 manuscript
addpath('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO')
mouse_num = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO\MouseNumInfoByRec.mat');
group = isVPA(mouse_num.mousenum);  
cd('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires');
folder_names = dir();
folder_names = arrayfun(@(n) folder_names(n).name,1:numel(folder_names),'UniformOutput',0);
temp = regexp(folder_names,'TrainingFit_Lambda4e-4_Kval','match');
temp = cellfun(@(x) ~isempty(x),temp,'UniformOutput',1);
folder_names(temp~=1)=[];
%loop through each K value and combile the data
num_motifs = NaN(sum(group==0),numel(folder_names));
num_k = NaN(1,numel(folder_names));
for cur_dir = 1:numel(folder_names)
   %Grab the k from the name   
   num_k(cur_dir) =str2double(folder_names{cur_dir}(regexp(folder_names{cur_dir},'Kval')+4:end));   
   fn = GrabFiles('TrainRepitoire_Block',0,{['Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\',folder_names{cur_dir}]});
   %remove saline
   [stats,~] = CompileStats(fn(group==0),{'ExpVarLoad'},[],'none');
   num_motifs(:,cur_dir) = arrayfun(@(n) numel(stats(n).ExpVarLoad),1:numel(stats),'UniformOutput',1); 
end
%sort everything by increasing k
[~,idx] = sort(num_k,'ascend');
num_k = num_k(idx);
num_motifs = num_motifs(:,idx);

%% PLot Figure
fp = fig_params;
figure; hold on; 

ci = bootci(100,@nanmedian,num_motifs);
avg_num = nanmedian(num_motifs);

shadedErrorBar(num_k,avg_num,cat(1,avg_num-ci(1,:),ci(2,:)-avg_num),'lineProps',{'color',fp.c_discovery,...
    'linewidth',fp.dl_line_width},'patchSaturation',fp.dl_alpha);
plot(num_k,avg_num,'o','color',fp.c_discovery,'markersize',6,'LineWidth',1.5)

%add line showing chosen length
line([28 28],[0 40],'linewidth',fp.dl_line_width,'linestyle','--','color',[0.1 0.1 0.1])
ylim([1,22])

%add titles and format
xlabel({'Maximum Number of';'Allowed Motifs'});
ylabel({'Number of Discovered';'Motifs'});
set(gca,'xtick',[1,15,30,45]);
xlim([num_k(1),num_k(end)]);

setFigureDefaults;
set(gca,'position',[3,5,6,6])
set(gcf,'position',[680   400  450   600]);

%%
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\ManuscriptRevisionFigures_currentbio';
handles = get(groot, 'Children');
saveCurFigs(handles,'-svg','K_sweep',savedir,1);
close all;

    
