function ExampleDataFigure(block)

%Camden MacDowell 2019 
if nargin <1
   block = [32];
end
savefigs = 1; 

%load the original data
orig = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DFF_Data\AllData_binned_SmallMask_3x_2minTraining.mat');

for cur = 1:numel(block)
    close all    
    cd('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28');
    savedir = ['C:\Users\macdo\OneDrive\Buschman Lab\AnalysisCode_Repository\Mesoscale Network Dynamics 2019 Analyses\Individual_Repitoire_Images_lesssmooth\' sprintf('Rec%d',block(cur))];
    if ~exist(savedir)
        mkdir(savedir);
    end
    
    data_orig = conditionDffMat(orig.data{block(cur)}',orig.nanpxs{block(cur)});
    [nX,nY,nT] = size(data_orig);
    nP =nX*nY;
    data_orig(isnan(data_orig))=0;
    data_orig = imgaussfilt3(data_orig,[1 1 0.1]);
    data_orig = reshape(data_orig,[nP,nT]);
    
    %Load the factorized data
    temp = load(sprintf('TrainRepitoire_block_%d.mat',block(cur)));
    
    %reconstruct data
    data_recon = helper.reconstruct(temp.w,temp.H);
    data_recon = conditionDffMat(data_recon',orig.nanpxs{block(cur)});
    data_recon(isnan(data_recon))=0;
    data_recon = imgaussfilt3(data_recon,[1 1 0.1]);
    data_recon = reshape(data_recon,[nP,nT]);

    %First plot the reconstruction
    figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
    imagesc(temp.Xhat)
    cmap = flipud(gray); 
    colormap(cmap)
    caxis([0 0.02]);    
    set(gca, 'ydir', 'reverse')
    axis off
    
    %Now the original
    figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
    data_orig_active = data_orig;
    data_orig_active(orig.nanpxs{block(cur)},:)=[];
    imagesc(data_orig_active)
    cmap = flipud(gray); 
    colormap(cmap)
    caxis([0 0.02]);    
    set(gca, 'ydir', 'reverse')
    axis off
    
    %Save off both as svgs
    if savefigs
       handles = get(groot, 'Children');
       saveCurFigs(handles,'-dsvg',sprintf('FullMatrix_block%d',block(cur)),savedir,0);
       close all
    end

    %Now plot the Hs on top of the reconstruction
    figure('units','normalized','Position',[0 0 1 1]);
    SimpleWHPlot(temp.w,temp.H);
	caxis([0 0.02])
    if savefigs
       handles = get(groot, 'Children');
       saveCurFigs(handles,'-dsvg',sprintf('SimpleWHPlotFull_block%d',block(cur)),savedir,0);
       close all
    end
    
    for i = 1:size(temp.W,2)
        figure();
        SimpleWHPlot(temp.W(:,i,:),temp.H(i,:));
        title(sprintf('Motif%d',i));
        if savefigs
           handles = get(groot, 'Children');
           saveCurFigs(handles,'-dsvg',sprintf('SimpleWHPlotFull_block%d_trace%d',block(cur),i),savedir,0);
           close all
        end
    end

    data = temp.W; 
    nL = size(data,3);
    %Now plot the individual motifs
    data_smooth = zeros(size(data));
    for cur_k = 1:size(data,2)
        temp = squeeze(data(:,cur_k,:));
        temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:))); %normalize
        temp = reshape(temp,[nX nY nL]);
        temp = imgaussfilt3(temp,[1, 1, 0.1]);   
        temp = reshape(temp,[nX*nY,nL]);
        data_smooth(:,cur_k,:) = temp; 
    end
    FlowAnalysis(data_smooth(:,[3,5,9],:),'SaveDir',savedir,'caxis',[0 99]);
    
    %Now do the same with static network
    Comparison = toggleModFlag(data_smooth,[],1);
    FlowAnalysis(Comparison(:,[3,5,9],:),'SaveDir',pwd,'caxis',[0 99],'TrigSTD',1000,'NameStr','Static');
    Comparison = toggleModFlag(data_smooth,[],2); %scale to motif intensity
    Resid = data_smooth-Comparison;
    FlowAnalysis_Residuals(Resid(:,[3,5,9],:),zeros(size(data_smooth,2)),'caxis',[0 99],'SaveDir',savedir);
    
    %code to plot a specific time series in frames
    savedir_recon = [savedir filesep 'FrameSubsetRecon1'];
    if ~exist(savedir_recon)
        mkdir(savedir_recon);
    end    
    FlowAnalysis(data_recon(:,1432:1452),'SaveDir',savedir_recon,'caxis',[0 99],'TrigSTD',100); %timepoint for example motif 1 = 1432:1452

    savedir_orig = [savedir filesep 'FrameSubsetOrig1'];
    if ~exist(savedir_orig)
        mkdir(savedir_orig);
    end    
    FlowAnalysis(data_orig(:,1432:1452),'SaveDir',savedir_orig,'caxis',[0 99],'TrigSTD',100); %timepoint for example motif 1 = 1432:1452
    
    savedir_recon = [savedir filesep 'FrameSubsetRecon2'];
    if ~exist(savedir_recon)
        mkdir(savedir_recon);
    end    
    FlowAnalysis(data_recon(:,1176:1187),'SaveDir',savedir_recon,'caxis',[0 99],'TrigSTD',100); %timepoint for example motif 1 = 1432:1452

    savedir_orig = [savedir filesep 'FrameSubsetOrig2'];
    if ~exist(savedir_orig)
        mkdir(savedir_orig);
    end    
    FlowAnalysis(data_orig(:,1176:1187),'SaveDir',savedir_orig,'caxis',[0 99],'TrigSTD',100); %timepoint for example motif 1 = 1432:1452
    
    
end

