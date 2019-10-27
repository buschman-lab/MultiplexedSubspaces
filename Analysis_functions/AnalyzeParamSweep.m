function [ev_norm, param] = AnalyzeParamSweep(in_fns,target_param)

%Outputs explained varainace and target param values 
%rows are runs, columns are param values within each run 
%currently target param is hardcoded, but this will change to just
%accepting field name 

%load data
temp = cellfun(@load, in_fns, 'UniformOutput', 0); 
temp = cellfun(@(x) x.paramsweep, temp, 'UniformOutput',0);

%get explained variance
fields = fieldnames(temp{1});
idx = cellfun(@(x) strcmp(x,'ExpVar_all'),fields,'UniformOutput',0);
idx = ([idx{:}]==1);
ev = cellfun(@(x) cat(2,x.(fields{idx})), temp,'UniformOutput',0);
ev_norm = cellfun(@(x) x/max(x), ev,'UniformOutput',0);
ev_norm = cat(1,ev_norm{:});

%get the param of interest
idx = cellfun(@(x) strcmp(x,target_param),fields,'UniformOutput',0);
idx = ([idx{:}]==1);
param = cellfun(@(x) cat(2,x.(fields{idx})), temp,'UniformOutput',0);
param = cat(1,param{:});


end %function 



















