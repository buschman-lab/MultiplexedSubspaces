function Plot_ExampleSTNMFofMotifs()
% Camden MacDowell - timeless

savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\ManuscriptRevisionFigures_currentbio\stnmfDecompositionOfMotifs';
if ~exist(savedir)
    mkdir(savedir)
end


%% on a per motif basis, plot the cummulative explained varaince of the unique number of spatial and temporal patterns
%If we are finding different ones for each motif, then usign the top 1
%dimension will capture most varaince. If we are mixing and matching them
%then I will take more. 
%sequentially add dimensions in A for each motif. Combine together and get
%explained variance. 

mouse_num = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO\MouseNumInfoByRec.mat');
group = isVPA(mouse_num.mousenum); 
n_dim = 40;    
rec_list = find(group==0);
pev = NaN(numel(rec_list),n_dim);
number_unique_across_motif= NaN(numel(rec_list),n_dim,2);
for cur_rec = 1:10%numel(rec_list)
    fprintf('\n\tWorking on recording %d of %d',cur_rec,numel(rec_list));
    temp = load(['Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\SpaceTimeNMF_Comparison\SpaceTimeNMF_xcorr_alligned_trackingfit_rankA\',...
        sprintf('block_%d.mat',rec_list(cur_rec))],'Ws','Wt','A');
    Ws = temp.Ws;
    Wt = temp.Wt;
    A = temp.A;
    %load the motifs
    temp = load(['Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\TrainRepitoire_block_',...
        num2str(rec_list(cur_rec)),'.mat'],'w');
    w= AlignMotifs(temp.w); 
    X = reshape(w,size(w,1),size(w,2)*size(w,3));
        
    for cur_dim = 1:n_dim     
        A_sparse = zeros(size(A));
%         number_unique = NaN(size(w,2),2);
        for cur_motif = 1:size(w,2)
            A_motif = squeeze(A(:,:,cur_motif));    
            %get the top x relationships in A
            idx = NaN(cur_dim,3);
            for i = 1:cur_dim
                [idx(i,3), ind] = max(A_motif(:));
                [idx(i,1),idx(i,2)] = ind2sub(size(A_motif),ind);
                A_motif(idx(i,1),idx(i,2))=0;
                %just keep the best relationships     
                A_sparse(idx(i,1),idx(i,2),cur_motif) = idx(i,3);
            end   
            if cur_motif == 1
               list_space_time = idx(:,1:2);                
            else
               list_space_time = cat(1,list_space_time,idx(:,1:2));
            end            
%             number_unique(cur_motif,:) = [numel(unique(idx(:,2))),numel(unique(idx(:,1)))];
        end %motif iteration
        
        

        %reconstruct the data
        X_recon = stNMF_Reconstruct(A_sparse,Wt,Ws);
        X_recon = X_recon';
        pev(cur_rec,cur_dim) = CalculateExplainedVariance(X,X-X_recon)*100;
        number_unique_across_motif(cur_rec,cur_dim,:) = [numel(unique(list_space_time(:,2))),numel(unique(list_space_time(:,1)))];
        

%         close; figure; hold on; 
%         subplot(2,1,1); imagesc(X_recon); ylabel('pixels'); xlabel('timepoints (concatenated motifs)');
%         subplot(2,1,2); imagesc(X);ylabel('pixels'); xlabel('timepoints (concatenated motifs)');
%         sgtitle(sprintf('PEV = %0.2g%%',pev(cur_rec,cur_dim)));
%         handles = get(groot, 'Children');
%         saveCurFigs(handles,'-svg',sprintf('ExampleTrace_dimensions_%d',cur_dim),savedir,1);        
        
    end %dim iteration
end %rec iteration

%% Plot the explained variance 

figure('position',[680  200  600   600]); hold on
pev_avg = nanmedian(pev);
ci = bootci(1000,@nanmedian,pev);
shadedErrorBar(1:numel(pev_avg),pev_avg,cat(1,pev_avg-ci(1,:),ci(2,:)-pev_avg),'lineprops',{'-','color',...
[0.25 0.25 0.25],'linewidth',2},'transparent',1,'patchSaturation',0.2);
ylim([0 90]);
ylabel({'Percent Explained Variance'})
xlabel({'Dimensions in A (per motif!)'})
setFigureDefaults;
set(gca,'position',[2,4,6,6])
handles = get(groot, 'Children');
%%
saveCurFigs(handles,'-svg','NumberOfDimensionsPerMotif',savedir,1);        

%% Plot the number of unique spatial and temporal patterns needed across each motif per epoch
figure('position',[680  200  600   600]); hold on
pev_avg = nanmedian(pev);
ci = bootci(1000,@nanmedian,pev);
shadedErrorBar(nanmedian(number_unique_across_motif(:,:,1)),pev_avg,cat(1,pev_avg-ci(1,:),ci(2,:)-pev_avg),'lineprops',{'-','color',...
[0.25 0.25 0.25],'linewidth',2},'transparent',1,'patchSaturation',0.2);
shadedErrorBar(nanmedian(number_unique_across_motif(:,:,2)),pev_avg,cat(1,pev_avg-ci(1,:),ci(2,:)-pev_avg),'lineprops',{'-','color',...
[0.25 0.25 0.25],'linewidth',2},'transparent',1,'patchSaturation',0.2);
ylim([0 90]);
ylabel({'Percent Explained Variance'})
xlabel({'Dimensions in A (per motif!)'})
setFigureDefaults;
set(gca,'position',[2,4,6,6])
handles = get(groot, 'Children');

%%%
%%

block = 32; 
motif = 9;
temp = load(['Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\SpaceTimeNMF_Comparison\SpaceTimeNMF_xcorr_alligned_trackingfit_rankA\',sprintf('block_%d.mat',block)],'Ws','Wt','A');
Ws = temp.Ws;
Wt = temp.Wt;
A = temp.A;

%load the motifs
temp = load(['Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\TrainRepitoire_block_' num2str(block),'.mat'],'w');
w= AlignMotifs(temp.w); 


A_motif = squeeze(A(:,:,motif));
nX = 20; %top x to plot
idx = NaN(nX,3);
for i = 1:nX
    [idx(i,3), ind] = max(A_motif(:));
    [idx(i,1),idx(i,2)] = ind2sub(size(A_motif),ind);
    A_motif(idx(i,1),idx(i,2))=0;
end

%% reconstruct each space and time relationship

for i = 1:nX
   A_temp = zeros(size(A));  
   A_temp(idx(i,1),idx(i,2),motif)=idx(i,3);
   
   X_recon = stNMF_Reconstruct(A_temp,Wt,Ws);
   X_recon = reshape(X_recon',size(w,1),size(w,2),size(w,3));
    
end
    
%%
X = reshape(w,size(w,1),size(w,2)*size(w,3));
X = reshape(X,size(w,1),size(w,2),size(w,3));
w_motif = squeeze(w(:,motif,:));


for i = 1:16
close; imagesc(squeeze(X(:,i,:)),[0 0.25]);
pause();
end


%%

































