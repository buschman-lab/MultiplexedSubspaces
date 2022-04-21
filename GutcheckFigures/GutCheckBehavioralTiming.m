function GutCheckBehavioralTiming(cur_rec,fast_load)
%Camden - timeless


%xcorr to make sure that parsing of behavioral timing is correct

if fast_load==0
    [~,~,~,EphysPath,~] = LoadDataDirectories(cur_rec);
    [st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe'); 
else %preprocessed
    [rec_name,~,~,EphysPath] = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDmotif1.mat'],0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA'}); %first motif has all that you need
    temp = cellfun(@(x) load(x,'st_mat','st_depth'),fn);  
    st_depth = temp.st_depth;
    st_mat = temp.st_mat; 
end

%normalize spiking rate, 
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);
st_norm = cellfun(@(x) x(1:2:end,:)+x(2:2:end,:),st_norm,'UniformOutput',0);

%average neural activity across everywhere
neu = nanmean(cat(2,st_norm{:}),2);

%load motor activity
[me,label] = AnalyzeBehavior(cur_rec,0);
me(1,:)=[];
me = (me-nanmean(me,2))';
[rho,lags] = xcorr(me(:,3), neu-nanmean(neu),30*(15/2),'normalized');
plot(lags,rho);

temp = circshift(me(:,3),8);
[rho,lags] = xcorr(temp, neu-nanmean(neu),30*(15/2),'normalized');
plot(lags,rho);


end


