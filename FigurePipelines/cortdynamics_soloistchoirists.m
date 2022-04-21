
%Camden MacDowell

%% Different neuron firing patterns be influenced by the normalization?
ImpactOfNormalization();

%% Load data
data = LoadSubspaceData('in_grouped');
datams = LoadSubspaceData('in_grouped_meansubtract');

%% compare the reduced rank regression performance across dimensions
[perf,area_label] = LoadVariable(data,'rel_performance',[],[],'mean');
perf = squeeze(nanmean(perf,1));
perfms = LoadVariable(datams,'rel_performance',[],[],'meansubtract');
perfms = squeeze(nanmean(perfms,1));

figure; hold on; 
col = fp.c_area; 
for cur_a = 1:8
    x = squeeze(perf(:,cur_a,:)-perfms(:,cur_a,:));
    shadedErrorBar(1:30,nanmean(x),sem(x,1),'lineprops',{'color',[col(cur_a,:),0.25],'linewidth',2});
end
ylabel('\Delta r^2')
xlabel('Subspace Dimension')
title({'Div vs subtract','Relative performance'},'fontweight','normal') 
fp.FormatAxes(gca); box on; grid on;
legend(area_label,'location','bestoutside')
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])

%% Are the beta weights between both groups correlated across dimensions? 
cur_a = 8; 
ndim = 10; 
cur_rec = 1; 
area_name = 'VIS';
x = NaN(14,cur_d);
for cur_d = 1:ndim
    beta = LoadVariable(data(cur_rec),'rrr_beta',area_name,cur_d,'mean');
    betams = LoadVariable(datams(cur_rec),'rrr_beta',area_name,cur_d,'meansubtract');
    %get correlation in betas per motif
    for cur_m = 1:14
        x(cur_m,cur_d) = corr(beta(cur_m,:)',betams(cur_m,:)','rows','complete');
    end
end
figure; hold on; 
shadedErrorBar(1:ndim,nanmean(x),sem(x,1),'lineprops',{'color',[col(cur_a,:),0.25],'linewidth',2});
ylabel('Rho')
xlabel('Subspace Dimension')
title({'Beta similarity'},'fontweight','normal') 
fp.FormatAxes(gca); box on; grid on;
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])


%% What are the properties of neurons driving these generalized dimensions? 
% plot the trial mean versus baseline (not residuals) for the mean div
cur_rec = 1; 
[tmean,area_label] = LoadVariable(datams(cur_rec),'trial_mean',area_name,[],'mean'); %%CAMDEN, Don't load it like this: won't work for the recs 3, and 4
tbase = LoadVariable(datams(cur_rec),'trial_baseline',area_name,[],'mean');
beta = LoadVariable(datams(cur_rec),'rrr_beta',area_name,5,'mean');
figure; hold on;
plot3(tmean(:),tbase(:),beta(:),'marker','o','linestyle','none')
xlabel('trial FR'); ylabel('baseline FR'); zlabel('beta weight')

%% 
cur_rec = 1; 
[tmean,area_label] = LoadVariable(data(cur_rec),'trial_mean',area_name,[],'mean');
tbase = LoadVariable(data(cur_rec),'trial_baseline',area_name,[],'mean');
beta = LoadVariable(data(cur_rec),'rrr_beta',area_name,5,'mean');
figure; hold on;
plot3(tmean(:),tbase(:),beta(:),'marker','o','linestyle','none')
xlabel('trial FR'); ylabel('baseline FR'); zlabel('beta weight')

%% get the spike triggered firing rate
%compute the spike-triggered relationship of each neuron relative 
spikes = LoadSubspaceData('in_grouped_raw');

area_label = [];
stPC = [];
rho = [];
for cur_rec = 1:6
    %parse activity per parent region 
    [area_val, area_label{cur_rec}] = ParseByArea(cat(2,spikes{cur_rec}.st_norm{:})',spikes{cur_rec}.neu_area,'general');

    %clean up areas %third input is the min # of spikes to keep area
    [area_val, area_label{cur_rec}] = CleanUpAreas(area_val, area_label{cur_rec}, 10); 

    %compute stPC population coupling
    [stPC{cur_rec}] = cellfun(@(x) PopulationCoupling(x,'stPC'), area_val,'UniformOutput',0);
    [rho{cur_rec}] = cellfun(@(x) PopulationCoupling(x,'rho'), area_val,'UniformOutput',0);
end 

%get MUA ID
ap_clusters = [];
for cur_rec = 1:6
    [~,~,~,EphysPath,~] = LoadDataDirectories(cur_rec);
    EphysPath = load(EphysPath);
    EphysPath = EphysPath.opts; 
    [ap_clusters{cur_rec},~] = CreateMasks(EphysPath); 
    [ap_clusters{cur_rec}, temp_area] = ParseByArea(cat(2,ap_clusters{cur_rec}{:})',spikes{cur_rec}.neu_area,'general');
    idx = ismember(temp_area,area_label{cur_rec});
    ap_clusters{cur_rec}(idx==0)=[];
end

[~,all_areas] = LoadVariable(data(1),'rrr_beta',[],5,'mean');
%%
%per motif and per recording, get the relationship between coupling and beta weight
rho = NaN(6,8,14,10);
for cur_rec = 1:6
    for cur_area = 1:numel(all_areas)
        area_name = all_areas{cur_area};
        idx = find(strcmp(area_label{cur_rec},area_name)==1);
        if ~isempty(idx)
            x = stPC{cur_rec}{idx};
            mask = ap_clusters{cur_rec}{idx};
            for cur_d = 1:10
                beta = LoadVariable(data,'rrr_V',area_name,cur_d,'mean'); 
                beta = squeeze(beta(cur_rec,:,:));
                beta(:,sum(isnan(beta),1)==14)=[];
%                 beta = abs(beta);
                for cur_m = 1:14
%                     y = beta(cur_m,mask==1)'; 
%                     xx = x(mask==1);    
                    y = beta(cur_m,:)'; 
                    xx = x(:);    
                    rho(cur_rec,cur_area,cur_m,cur_d) = corr(y,xx,'type','pearson');
                    %also save off MUA vs NUE betas 
                end
            end
        end
    end
end

%% Plot the correlation across areas and dimensions
%plot per area
figure; hold on; 
col = fp.c_area;
for cur_a = 1:8
    temp = fisherZ((rho));
%     temp = squeeze(nanmean(temp,1));
%     temp = squeeze(temp(cur_a,:,:));
    temp = squeeze(temp(:,cur_a,:,:));
    temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
    shadedErrorBar(1:10,nanmean(temp),sem(temp,1),'lineprops',{'color',[col(cur_a,:),0.25],'linewidth',2});
end
legend(all_areas)


%to do | 
%plot the above for each area opposite of the share dimensionality for that
%area. 
%plot correlation between generalization of each area and dimension versus
%the soloist/choir
%same as above for MUAs
%put together all as a figure 3 and provide to Tim

%% Plot local coupling versus beta weights
%per motif and per recording, get the relationship between coupling and beta weight
rho = NaN(6,8,14,10);
for cur_rec = 1:6
    for cur_area = 1:numel(all_areas)
        area_name = all_areas{cur_area};
        idx = find(strcmp(area_label{cur_rec},area_name)==0);
        if ~isempty(find(strcmp(area_label{cur_rec},area_name)==1))
            x = cat(1,stPC{cur_rec}{idx});
            mask = cat(1,ap_clusters{cur_rec}{idx});
            for cur_d = 1:10
                beta = LoadVariable(data,'rrr_beta',area_name,cur_d,'mean'); 
                beta = squeeze(beta(cur_rec,:,:));
                beta(:,sum(isnan(beta),1)==14)=[];
%                 beta = abs(beta);
                for cur_m = 1:14
                    y = beta(cur_m,mask==1)'; 
                    xx = x(mask==1);    
                    rho(cur_rec,cur_area,cur_m,cur_d) = corr(y,xx,'type','pearson');
                    %also save off MUA vs NUE betas 
                end
            end
        end
    end
end

%% Plot the correlation across areas and dimensions
%plot per area
figure; hold on; 
col = fp.c_area;
for cur_a = 1:8
    temp = fisherZ((rho));
%     temp = squeeze(nanmean(temp,1));
%     temp = squeeze(temp(cur_a,:,:));
    temp = squeeze(temp(:,cur_a,:,:));
    temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
    shadedErrorBar(1:10,nanmean(temp),sem(temp,1),'lineprops',{'color',[col(cur_a,:),0.25],'linewidth',2});
end
legend(all_areas)



%% Then show that the neurons being weighted most differently between the two conditions are the ones that fit this criteria



















