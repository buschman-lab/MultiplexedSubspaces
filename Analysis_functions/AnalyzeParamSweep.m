function [ev_norm, params] = AnalyzeParamSweep(in_fns,target_params)

%Outputs explained varainace and target param values 
%rows are runs, columns are param values within each run 
%currently target param is hardcoded, but this will change to just
%accepting field name 

%load data
data = cellfun(@load, in_fns, 'UniformOutput', 0); 
data = cellfun(@(x) x.paramsweep, data, 'UniformOutput',0);

%get explained variance
fields = fieldnames(data{1});
idx = cellfun(@(x) strcmp(x,'ExpVar_all'),fields,'UniformOutput',0);
idx = ([idx{:}]==1);
ev = cellfun(@(x) cat(2,x.(fields{idx})), data,'UniformOutput',0);
ev_norm = cellfun(@(x) x/max(x), ev,'UniformOutput',0);
ev_norm = cat(1,ev_norm{:});

%get the param of interest
params = NaN(numel(target_params),size(ev_norm,2));
for i = 1:numel(target_params)
    idx = cellfun(@(x) strcmp(x,target_params{i}),fields,'UniformOutput',0);
    idx = ([idx{:}]==1);
    temp = cellfun(@(x) cat(2,x.(fields{idx})), data,'UniformOutput',0);
    temp = cat(1,temp{:});
    if iscell(temp); temp = cell2mat(temp); end
    params(i,:) = nanmean(temp); %get the average parameter across all trials
end

end %function 



















