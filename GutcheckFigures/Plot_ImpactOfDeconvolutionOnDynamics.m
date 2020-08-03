function Plot_ImpactOfDeconvolutionOnDynamics()

data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep';

file_list = GrabFiles('.mat',0,{data_dir});

%load data
temp = cellfun(@(x) load(x,'W','data_train','stats_test','stats_train'),file_list,'UniformOutput',0);
W = cellfun(@(x) x.('W'),temp,'UniformOutput',0);
stats_test = cellfun(@(x) x.('stats_test'),temp,'UniformOutput',0);
stats_train = cellfun(@(x) x.('stats_train'),temp,'UniformOutput',0);
data_train = temp{1}.data_train;

%get the smoothing parameter from the file name:
smoothval = cellfun(@(x) str2num(x((regexp(x,'smooth','end')+1):regexp(x,'.mat','start'))),file_list,'UniformOutput',1);
        
fp = fig_params_vpa;

dynamicness = zeros(1,numel(W));
stvar = zeros(1,numel(W));
for i = 1:numel(W)
    temp = W{i};
    dynamicness(i) =  1 - fisherInverse(nanmean(arrayfun(@(n) nanmean(fisherZ(corr(squeeze(temp(:,n,:)),nanmean(squeeze(temp(:,n,:)),2)))),1:size(temp,2),'UniformOutput',1)));   
    stvar(i) = nanmean(arrayfun(@(n) nanvar(reshape(squeeze(temp(:,n,:)),1,numel(temp(:,n,:)))),1:size(temp,2),'UniformOutput',1));
end

train_pev = cellfun(@(x) x.pev,stats_train,'UniformOutput',1);
test_pev = cellfun(@(x) x.pev,stats_test,'UniformOutput',1);
num_motifs = cellfun(@(x) x.n_motifs,stats_train,'UniformOutput',1);
lambda = cellfun(@(x) x.lambda,stats_train,'UniformOutput',1);
    
%%plot figures;
xval = unique(smoothval);

%dynamicsness
figure; hold on; 
y = arrayfun(@(x) nanmean(dynamicness(smoothval==x)),xval,'UniformOutput',1);
plot(xval,y,'marker','o','linewidth',2,'color','k');
xlabel('deconvolution smoothing')
ylabel('dissimilarity to mean')
title('impact on dynamicsness')

%variance
figure; hold on; 
y = arrayfun(@(x) nanmean(stvar(smoothval==x)),xval,'UniformOutput',1);
plot(xval,y,'marker','o','linewidth',2,'color','k');
xlabel('deconvolution smoothing')
ylabel('variance')
title('impact on spatiotemporal variance')

%train pev
figure; hold on; 
y = arrayfun(@(x) nanmean(train_pev(smoothval==x)),xval,'UniformOutput',1);
plot(xval,y,'marker','o','linewidth',2,'color','k');
xlabel('deconvolution smoothing')
ylabel('train pev')
title('impact on quality of fit')

%test pev
figure; hold on; 
y = arrayfun(@(x) nanmean(test_pev(smoothval==x)),xval,'UniformOutput',1);
plot(xval,y,'marker','o','linewidth',2,'color','k');
xlabel('deconvolution smoothing')
ylabel('test pev')
title('impact on generalizability of fit')

%num_motifs
figure; hold on; 
y = arrayfun(@(x) nanmean(num_motifs(smoothval==x)),xval,'UniformOutput',1);
plot(xval,y,'marker','o','linewidth',2,'color','k');
xlabel('deconvolution smoothing')
ylabel('num motifs')
title('impact on number of motifs')

%lambda 
figure; hold on; 
y = arrayfun(@(x) nanmean(lambda(smoothval==x)),xval,'UniformOutput',1);
plot(xval,y,'marker','o','linewidth',2,'color','k');
xlabel('deconvolution smoothing')
ylabel('autofit lambda')
title('impact on overlap penalty')

%% Now just visualize the discovered motifs
%load one example animal
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep';
file_list = GrabFiles('Mouse9069\w*',0,{data_dir});
%Grabfiles is being weird so just keep gif ones
temp = cellfun(@(x) ~isempty(regexp(x,'gif','ONCE')),file_list,'UniformOutput',1);
file_list(temp)=[];

[path,name] = cellfun(@(x) fileparts(x),file_list,'UniformOutput',0);
%load data
temp = cellfun(@(x) load(x,'W','nanpxs'),file_list,'UniformOutput',0);
W = cellfun(@(x) x.('W'),temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.('nanpxs'),temp,'UniformOutput',0);
for i = 1:numel(W)
    temp = W{i};
    for j = 1:size(temp,2)
        motif = conditionDffMat(squeeze(temp(:,j,:))',nanpxs{i},[],[68 68 size(temp,3)]);
        MotifToGif(motif,[path{i},filesep, name{i}, sprintf('_motif%d.gif',j)])
    end
end

%% Compare with 2std threshold
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep\Comparisons using 2std threshold';
file_list = GrabFiles('Mouse9069\w*',0,{data_dir});
%Grabfiles is being weird so just keep gif ones
temp = cellfun(@(x) ~isempty(regexp(x,'gif','ONCE')),file_list,'UniformOutput',1);
file_list(temp)=[];

[path,name] = cellfun(@(x) fileparts(x),file_list,'UniformOutput',0);
%load data
temp = cellfun(@(x) load(x,'W','nanpxs'),file_list,'UniformOutput',0);
W = cellfun(@(x) x.('W'),temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.('nanpxs'),temp,'UniformOutput',0);
for i = 1:numel(W)
    temp = W{i};
    for j = 1:size(temp,2)
        motif = conditionDffMat(squeeze(temp(:,j,:))',nanpxs{i},[],[68 68 size(temp,3)]);
        MotifToGif(motif,[path{i},filesep, name{i}, sprintf('_motif%d.gif',j)])
    end
end
%% Compare with 1std threshold
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep\Comparisons using 1std threshold';
file_list = GrabFiles('Mouse9069\w*',0,{data_dir});
%Grabfiles is being weird so just keep gif ones
temp = cellfun(@(x) ~isempty(regexp(x,'gif','ONCE')),file_list,'UniformOutput',1);
file_list(temp)=[];

[path,name] = cellfun(@(x) fileparts(x),file_list,'UniformOutput',0);
%load data
temp = cellfun(@(x) load(x,'W','nanpxs'),file_list,'UniformOutput',0);
W = cellfun(@(x) x.('W'),temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.('nanpxs'),temp,'UniformOutput',0);
for i = 1:numel(W)
    temp = W{i};
    for j = 1:size(temp,2)
        motif = conditionDffMat(squeeze(temp(:,j,:))',nanpxs{i},[],[68 68 size(temp,3)]);
        MotifToGif(motif,[path{i},filesep, name{i}, sprintf('_motif%d.gif',j)])
    end
end



%% Plot an example clip from the training data using each deconvolution, 1std, and 2std

%1std
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep\Comparisons using 1std threshold';
file_list = GrabFiles('Mouse9069\w*',0,{data_dir});
%Grabfiles is being weird so just keep gif ones
temp = cellfun(@(x) ~isempty(regexp(x,'gif','ONCE')),file_list,'UniformOutput',1);
file_list(temp)=[];

[path,name] = cellfun(@(x) fileparts(x),file_list,'UniformOutput',0);
%load data
temp = cellfun(@(x) load(x,'data_test','nanpxs'),file_list,'UniformOutput',0);
data_test = cellfun(@(x) x.('data_test'),temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.('nanpxs'),temp,'UniformOutput',0);
for i = 1:numel(data_test)   
    clip = conditionDffMat(data_test{i}',nanpxs{i},[],[68 68 size(data_test{i},2)]);
    clip = clip(:,:,1:500);
    MotifToGif(clip,[path{i},filesep, name{i}, '_1stdthreshold.gif'],'limit',99)
end

%2std
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep\Comparisons using 2std threshold';
file_list = GrabFiles('Mouse9069\w*',0,{data_dir});
%Grabfiles is being weird so just keep gif ones
temp = cellfun(@(x) ~isempty(regexp(x,'gif','ONCE')),file_list,'UniformOutput',1);
file_list(temp)=[];

[path,name] = cellfun(@(x) fileparts(x),file_list,'UniformOutput',0);
%load data
temp = cellfun(@(x) load(x,'data_test','nanpxs'),file_list,'UniformOutput',0);
data_test = cellfun(@(x) x.('data_test'),temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.('nanpxs'),temp,'UniformOutput',0);
for i = 1:numel(data_test)   
    clip = conditionDffMat(data_test{i}',nanpxs{i},[],[68 68 size(data_test{i},2)]);
    clip = clip(:,:,1:500);
    MotifToGif(clip,[path{i},filesep, name{i}, '_2stdthreshold.gif'],'limit',100)
end

%deconvolution
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep';
file_list = GrabFiles('Mouse9069\w*',0,{data_dir});
%Grabfiles is being weird so just keep gif ones
temp = cellfun(@(x) ~isempty(regexp(x,'gif','ONCE')),file_list,'UniformOutput',1);
file_list(temp)=[];

[path,name] = cellfun(@(x) fileparts(x),file_list,'UniformOutput',0);
%load data
temp = cellfun(@(x) load(x,'data_test','nanpxs'),file_list,'UniformOutput',0);
data_test = cellfun(@(x) x.('data_test'),temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.('nanpxs'),temp,'UniformOutput',0);
for i = 1:numel(data_test)   
    clip = conditionDffMat(data_test{i}',nanpxs{i},[],[68 68 size(data_test{i},2)]);
    clip = clip(:,:,1:500);
    MotifToGif(clip,[path{i},filesep, name{i}, '_deconv.gif'])
end

end
