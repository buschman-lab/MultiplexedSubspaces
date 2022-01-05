function CombinedEphysImagingVideo(ImgPath,BehavPath,EphysPath,ImgProbeLoc,EyeCrop,figtitle,savedir,behavstr,preprocessedFLAG)
%Camden MacDowell - timeless
%notes for future camden. The depth of probe insertions were reasonably
%consistent between recordings. Can tile the size of of the figure by that,
%which helps account for differences insertion depth. 

%% to do
%figure out the depth issue
% savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\RawDataVideos';
% figtitle = 'Mouse 332 Recording 1';
% ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined_processed.mat';
% ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_probe_coords.mat';
% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
% BehavPath = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Cam_1_20210607-144148.avi';
% % BehavPath2 = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Cam_0_20210607-144148.avi';
% preprocessedFLAG = 1; %if plotting the deconvolved or nondeconvolved data
% behavstr = [579,34048];
%Rectangle on the behavioral video to crop to visualize eyeball
% EyeCrop = [260,180,65,60];

%manual inspection revealed that the behavioral video is no sampling
%specifically at 60 hz. So take the start and stop time and space indices
behavidx = round(linspace(behavstr(1),behavstr(2),81000));
%adjust for the deconvolution
behavidx = behavidx(15+1:end-15);

%load imaging 
if preprocessedFLAG == 1
    load(ImgPath,'data_norm','nanpxs');  
    cvals = [0 8];
else
    load(ImgPath,'dff','nanpxs');
    data_norm = dff(1:2:end,:) + dff(2:2:end,:); clear dff;
    data_norm = data_norm'; clear dff;     
    data_norm = data_norm(:,15+1:end-15);
    cvals = [-8 8];
end
data = conditionDffMat(data_norm',nanpxs); 
probe_loc = load(ImgProbeLoc); probe_loc = probe_loc.probe_coords;

%load ephys
[~,opts,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe');
[st_mat,opts,st_vert_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','vert');
TileSize = [4,4,4,3];
%standard deviation normalized firing rate
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);
%reorganize from top to bottom (m,s,par,v)
p_order = [4,2,1,3]; %for colors
st_norm = st_norm(p_order);
st_depth = st_depth(p_order);
st_vert_depth = st_vert_depth(p_order);
win = 30; %window length in frames
prewin = 30; %window length in frames
fr_max = 5; %max (normalized) firing rate for scaling

%get the anatomical locations
ccf_path = fileparts(EphysPath);
ccf_path = load([ccf_path filesep 'probe_ccf_xvalidated.mat']);
ccf_path.probe_ccf = ccf_path.probe_ccf(p_order);
neu_area = arrayfun(@(n) MapAnatomicalLocation(ccf_path.st,ccf_path.probe_ccf(n),st_depth{n},1),1:numel(st_depth),'UniformOutput',0);
%inverse labels to match the depth plotting (inversed)
neu_area = cellfun(@(x)  x(linspace(numel(x),1,numel(x))), neu_area,'UniformOutput',0);
%switch to acroynm
neu_area = cellfun(@(x) AreaAcryonym(x,ccf_path.st), neu_area,'UniformOutput',0);
anotation_type = 'detail'; %detailed structures or parent structure names

%get the deconvolved roi and average spiking per probe
x = sqrt(size(data_norm,1)+numel(nanpxs)); 
r = 2;
offset = {[2,0],[0,-1],[1,-1],[0 0]};
%parse the radius around the probes
dff_probe = NaN(size(data,3),numel(probe_loc));
st_probe = NaN(size(data,3),numel(probe_loc));
for cur_probe = 1:numel(probe_loc)
    %create mask of radius around probe tip
    temp = probe_loc{cur_probe};
    mask = zeros(x,x);  
    temp(1,1)=temp(1,1)+offset{cur_probe}(1);
    temp(1,2)=temp(1,2)+offset{cur_probe}(2);
    mask(temp(1,2)-r:temp(1,2)+r,temp(1,1)-r:temp(1,1)+r)=1;
    %flatten
    mask = mask(:);
    mask(nanpxs)=[];        
    
    dff_probe(:,cur_probe) = nanmean(data_norm(mask==1,:),1); 
    st_probe(:,cur_probe) = nanmean(st_mat{cur_probe},2)/std(nanmean(st_mat{cur_probe},2));
end    
st_probe = st_probe(:,p_order);
dff_probe = dff_probe(:,p_order);

fp = fig_params_deconvolutionpaper;

%behavioral video
v = VideoReader(BehavPath);
if numel(EyeCrop)==5 %some are flipped, so add a throwaway 5th term to be passed into function
    flipimg = 1;
    EyeCrop = EyeCrop(1:4);
else
    flipimg=0;
end

%% Video 1 | Full data
col = getColorPalet(4);
% create the video writer with 15 fps
writerObj = VideoWriter([savedir filesep strrep(figtitle,' ','_') sprintf('_imgephysbehav%d.avi',preprocessedFLAG)]);
writerObj.FrameRate = 15;
open(writerObj);

try
    for cur_f = prewin+1:20000   
        close all; 
        figure('units','normalized','position',[0.2719 0.0657 0.6979 0.8398]); hold on; 
        set(gcf,'color','w')
        t=tiledlayout(4,6,'Padding','normal','TileSpacing','normal');
        title(t,figtitle)

        %plot the imaging frames
        nexttile([2,2]); hold on;
        temp = SpatialGaussian((data(:,:,cur_f)),3,1);
        imagesc(temp,cvals);
        set(gca,'ydir','reverse'); colormap(gca,'magma'); 
        ylim([6,62]); xlim([4,63]); axis equal    
        axis off            
        arrayfun(@(n) plot(probe_loc{n}(1,1)+offset{n}(1),probe_loc{n}(1,2)+offset{n}(2),'linewidth',1,'color',col(n,:),'marker','x','markersize',10),p_order,'UniformOutput',0)            
        
        rng(4);
        for cur_p = 1:4            
            c = distinguishable_colors(80); 
            nexttile([TileSize(cur_p),1]); hold on;
            imagesc(st_norm{cur_p}(cur_f-prewin:cur_f+win,:)',[0 fr_max]);  
            set(gca,'ydir','reverse')
            PlotProbeAnatomy(gca, neu_area{cur_p}, 0, anotation_type,1,2,c(randperm(size(c,1),size(c,1)),:));  
    %         colormap(flipud(gray))
            colormap(gca,flipud(gray));            
            yvals = floor(linspace(1,size(st_norm{cur_p},2),10));
%             set(gca,'YTick',yvals,'YTickLabel',st_depth{cur_p}(yvals),'YTickLabelRotation',45);               
            set(gca,'YTick','')
            if cur_p == 4
               c = colorbar('southoutside');
               ylabel(c,'Normalized Firing Rate')
            end
            %plot the current imaging frame line
            ylim([0,yvals(end)])
            ymax = get(gca,'ylim');
            plot([prewin+1,prewin+1],ymax,'linewidth',1.5,'color','r','linestyle','--'); %plus 1 so that it is on the current frame
            %plot a line at 600uM
%             sup_cort = find(st_vert_depth{cur_p}>=600,1,'first');
%             plot([0,win+prewin+1],[sup_cort,sup_cort],'linewidth',1.5,'color','b','linestyle','-'); %plus 1 so that it is on the current frame
            title(sprintf('probe %d',cur_p),'fontweight','normal','color',col(p_order(cur_p),:))
            xlim([0,win+prewin+1]);  
            xvals =get(gca,'xlim');
            set(gca,'xtick',[xvals(1) prewin+1 xvals(2)],'xticklabels',[-1*win/15,0 win/15]);            
            xlim([-10,win+prewin+1]);  
            set(gca,'XColor','w','YColor','w')
            fp.FormatAxes(gca)            
            box off
        end %probe loop    

        %behavioral video
        nexttile([2,2]); hold on;         
        if flipimg==1
            frame = flipud(read(v,behavidx(cur_f)));            
        else
            frame = read(v,behavidx(cur_f));
        end
        frame = imcrop(frame,[50,0,499,500]);
        eyeball = imresize(imadjust(im2gray(imcrop(frame,EyeCrop))),3);
        frame([(size(frame,1)+1-size(eyeball,1)):size(frame,1)],[1:size(eyeball,2)],:)=repmat(eyeball,[1,1,3]);
        imagesc(frame);  colormap(gca,'gray');
        set(gca,'ydir','normal'); 
        ylim([0 size(frame,1)]); xlim([0 size(frame,2)]); axis off
        axis equal    

        writeVideo(writerObj, getframe(gcf));
    end
catch
    %make sure the video closes
    close(writerObj);
end

% create the video writer with 15 fps
close(writerObj);


%% Video of the average traces
% create the video writer with 15 fps
writerObj = VideoWriter([savedir filesep strrep(figtitle,' ','_') sprintf('_avgfr%d.avi',preprocessedFLAG)]);
writerObj.FrameRate = 15;
open(writerObj);
for cur_f = prewin+1:20000   
    close all;
    figure('units','normalized','position',[0.2719 0.0657 0.6979 0.2148]); hold on; 
    set(gcf,'color','w')
    t=tiledlayout(1,5,'TileSpacing','compact');
    title(t,figtitle);
    %plot the imaging frames
    nexttile([1,1]); hold on;
    temp = SpatialGaussian((data(:,:,cur_f)),3,1);
    imagesc(temp,cvals);
    set(gca,'ydir','reverse'); colormap(gca,'magma'); 
    ylim([6,62]); xlim([6,65]); axis square
    axis off
    arrayfun(@(n) plot(probe_loc{n}(1,1),probe_loc{n}(1,2),'linewidth',1,'color',col(n,:),'marker','x','markersize',10),p_order,'UniformOutput',0)
    
    %plot the average activity from the probes and cranio fluorescence
    for cur_p = 1:4
        nexttile([1,1]); hold on;
        plot(dff_probe(cur_f-prewin:cur_f+win,cur_p),'k','linewidth',2,'linestyle','-');
        plot(st_probe(cur_f-prewin:cur_f+win,cur_p),'k','linewidth',2,'linestyle',':');
        title(sprintf('probe %d',cur_p),'fontweight','normal','color',col(p_order(cur_p),:))
        ylim([0 max(cat(1,st_probe(:),dff_probe(:)))]);
        plot([prewin+1,prewin+1],[0 max(cat(1,st_probe(:),dff_probe(:)))],'linewidth',1.5,'color','r','linestyle','--'); %plus 1 so that it is on the current frame        
        xlim([0,win+prewin+1]);
    end   
    writeVideo(writerObj, getframe(gcf));    
end

% create the video writer with 15 fps
close(writerObj);
close all;

end %function

