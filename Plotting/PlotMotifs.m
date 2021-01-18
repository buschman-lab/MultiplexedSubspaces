function PlotMotifs(data,varargin)
%Camden MacDowell 2020
% Plots motifs in a purty way. The figure is pretty. The code is not. 
%
% @INPUTS
% @data (required): an N x F x T stack. Will make T figures for each F. 
% individual image size is square: sqrt(N) x sqrt(N);. If data is N x T,
% the code will add an singleton F dimension. 
%
% opts (optional): see below
%
% If you want to run in a loop of folders, do:
% folders = dir('TrainingFit*');
% start = pwd;
% for cur = 1:numel(folders)
%     cd(folders(cur).name);
%     temp = load('Motifs.mat');
%     if exist([pwd filesep 'PlotMotifs'])==0
%         mkdir('PlotMotifs');
%     end
%     cd('PlotMotifs');
%     PlotMotifs(temp.W_clust_smooth);
% end
% 

close all

%Set options
%Define options
opts.SaveDir = pwd; %Save location
opts.Verbose = 1; %how chatty should we be?
opts.FontSize = 16;
opts.FontWeight = 'bold';
opts.NameStr = '';
opts.BregmaX = 1.97; %conversion ratio for finding location to plot bregma
opts.BregmaY = 2.27; %conversion ratio for finding locatio to plot bremga
opts.MaskDir = 'Z:\Rodent Data\Wide Field Microscopy\VPA_Mesoscale_Analysis\FigureMask.mat'; %the mask used for figure plotting
opts.caxis = [0 95]; %percentile of maximum intensity
opts.kernel = [];
opts.includeflow = 0; %add flow arrows?
opts.Percentile = 50; %percent of pixels to perform analysis on. 
opts.QuiverScale = 0.75; %arrow scaling factor for quiver
opts.FlowDownsample = 2; 
%Process optional inputs
opts = ParseOptionalInputs(opts,varargin);

%load mask;
mask = load(opts.MaskDir);
mask = mask.mask;

%break mask into l and r hemisphere
temp = mask;
temp(:,34:end) = 0;
hemi_mask{1} = temp; 
temp = mask;
temp(:,1:34) = 0;
hemi_mask{2} = temp; 

for cur_F = 1:size(data,2)
    cur_data = squeeze(data(:,cur_F,:));
    nX = sqrt(size(cur_data,1));
    nY = sqrt(size(cur_data,1));
    
    %Bregma coordinates
    bX = nX/opts.BregmaX;
    bY = nY/opts.BregmaY;
    
    cur_data = reshape(cur_data,nX,nY,size(cur_data,2));
    if ~isempty(opts.kernel)
        for i = 1:size(cur_data,3)
            cur_data(:,:,i) =  imgaussfilt(cur_data(:,:,i),opts.kernel,'filterdomain','spatial','FilterSize',5);
        end
    end
    cur_data(:,:,end+1) = zeros(nX,nY);
    

    %%
    climits = [prctile(cur_data(:),opts.caxis(1)),prctile(cur_data(:),opts.caxis(2))];
    
    %Loop through the individual images of the rec
    for cur_T = 1:size(cur_data,3)-1        
        %Create figure
        fig = figure();  hold on
        camroll(180)        
        
        %loop through each hemisphere
        for hemi = 1:2
            %Black background (to smooth the pixelated edges)
            cc = bwconncomp(hemi_mask{hemi},8);
            s = regionprops(cc,'Area','Centroid','ConvexHull');              
            fill(s(1).ConvexHull(:,1),s(1).ConvexHull(:,2),'k','LineWidth',2);
               
            cur_img = cur_data(:,:,cur_T);
            nxt_img = cur_data(:,:,cur_T+1);
            cur_img(~hemi_mask{hemi})=0;
            nxt_img(~hemi_mask{hemi})=0;
            
                        
            %Plot the raw data  
            fH = imagesc(cur_img); hold on

            %Make background transparent
            set(fH, 'AlphaData',hemi_mask{hemi});            
            
            %Draw a smooth border to smooth the pixelations
            plot(s(1).ConvexHull(:,1),s(1).ConvexHull(:,2),'Color','w','LineWidth',3);
    
            %Set figure parameters
            colormap magma
            axis square
            caxis(climits)
            axis off
            title(sprintf('Frame %d',cur_T),'Color','k','VerticalAlignment','baseline','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
                        %If suprathreshold intensity, do flow analysis
            if opts.includeflow
                %Flow analysis using HS method
                [u, v] = HS_flow(cur_img,nxt_img, 0.35, 2000);

                %Remove thresh_masked pxls
                thresh_mask = cur_img;
                thresh_mask(thresh_mask==0)=NaN;
                thresh_mask(thresh_mask<= prctile(thresh_mask,opts.Percentile,'all'))=0;
                thresh_mask(isnan(thresh_mask))=0;
                thresh_mask(thresh_mask~=0)=1;     
                u = u.*thresh_mask;
                v = v.*thresh_mask; 

                %Downsample and plot
                [gridX, gridY, u] = mean_downsample(u,opts.FlowDownsample);
                [~,~,v] = mean_downsample(v,opts.FlowDownsample);
                idx=(1:3:numel(gridX));%subsample (e.g. so you don't swamp with arrows)
                fQ = quiver(gridX(idx),gridY(idx),u(idx,idx),v(idx,idx),opts.QuiverScale,'color','g','LineWidth',1.7);          
            end %if trigger threshold 

        end %Hemi loop    
        %Addbregma
        scatter(bX,bY,600,'.','MarkerFaceColor',[1 0.2 .2],'MarkerEdgeColor',[1 0.2 .2]); 
        ylim([0 68]);
        xlim([0 68]);
        
        %Set Axes
        ylim([0 68]);
        xlim([0 68]);
        drawnow      
        
        set(gca,'FontSize',opts.FontSize,'FontWeight',opts.FontWeight);
        
    end %cur_T loop 
    
%%
    %Save off all files from the cur_F
    handles = get(groot, 'Children');
    handles = handles(sort(double(findall(0, 'type', 'figure'))));    
    filename = sprintf('%sMotifs_Group%d.png',opts.NameStr,cur_F);
    saveCurFigs(handles,'-png',filename,opts.SaveDir,1)
    close all    

% %     %Make Gif
%     MakeRepitoireGifs(data(:,cur_F,:),'Individual',1,'TrackPercentile',opts.TrackPercentile,...
%         'Track',1,'fnbase', sprintf('Group%d.gif',cur_F),'savedir',opts.SaveDir,...
%         'Combined',0,'TrigSTD',opts.TrigSTD);       
%     close all

end %cur_F loop

end %function end



