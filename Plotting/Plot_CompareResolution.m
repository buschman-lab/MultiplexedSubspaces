function fh = Plot_CompareResolution(in_fns,fh)

if nargin <2; fh = figure; end
set(0, 'currentfigure', fh);
hold on;

%load figure class
fp = fig_params; 

%set the x axis labels
labels = {};
target = 'numFactors'; %the target parameter to plot from the structure

%Outputs explained varainace and target param values 

%load data
data = cellfun(@load, in_fns, 'UniformOutput', 0); 
data = cellfun(@(x) x.stat, data, 'UniformOutput',0);

%parse data structure
fields = fieldnames(data{1});
idx = cellfun(@(x) strcmp(x,target),fields,'UniformOutput',0);
idx = ([idx{:}]==1);
param = cellfun(@(x) cat(2,x.(fields{idx})), data,'UniformOutput',0);
param = cat(1,param{:});
if iscell(param); param = cell2mat(param); end

%transpose
param = param';


%Convert to using gramm
%violin plot if enough points
if size(param,2)>20
    vp = CompareViolins(param,fp);
else %otherwise bar and jitter
    boxplot(param','color','k');
    set(findobj(gca,'type','line'),'linew',1.5)
    set(gca,'ylim',[0 42])
end
[p, tbl, ~] = anova1(param',[],'off');
AddSig(1,p,[1,3,41,41],2,2,1)
set(gca,'Xticklabel',{'136x136','68x68p','68x68 smooth'},'xticklabelrotation',45);
xlabel('Spatial Resolution');
ylabel('# Discovered Motifs');
setFigureDefaults;
set(gca,'position',[3,5,6,6])
title({'Number of Discovered Motifs';'is Robust to Spatial Resolution'},...
    'FontName','Arial','FontSize',16,'FontWeight','normal','position',[2 45]);
set(fh,'position',[680   400  450   600]);



%%
if savefigs
   handles = get(groot, 'Children');
   saveCurFigs(handles,'-svg','Within_between_expvar',savedir,1);
   close all
end




end %function 
















