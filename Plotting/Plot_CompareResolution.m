function fh = Plot_CompareResolution(in_fns,fh)

if nargin <2; fh = figure; end
set(0, 'currentfigure', fh);
hold on;

%load figure class
fp = fig_params; 

%set the x axis labels
labels = {}:
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
    shadedErrorBar(1:size(param,2),mean(param,1),sem(param,1),'lineProps',{'color',fp.c_discovery,...
        'linewidth',fp.dl_line_width},'patchSaturation',fp.dl_alpha);
else %otherwise bar and jitter
    


%add titles and format
xlabel('Motif Length (ms)');
ylabel({'Percent Explained Variance','(normalized)'});
fp.SetTitle(gca,{'Post-Hoc Validation 1s Motif Duration'});

fp.FormatAxes(gca)




end %function 



















