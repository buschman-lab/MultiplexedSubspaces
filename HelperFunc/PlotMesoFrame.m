function fH=PlotMesoFrame(frame,varargin)
%camden macdowell - timeless
%frame is a 68 x 68 pixel frame
%Set options
%Define options
opts.BregmaX = 1.97; %conversion ratio for finding location to plot bregma
opts.BregmaY = 2.27; %conversion ratio for finding locatio to plot bremga
opts.MaskDir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis\Preprocessing\brainoutline_small.mat'; %the mask used for figure plotting
opts.caxis = [0 95]; %percentile of maximum intensity
opts.caxis_flag = 0; %flag to make caxis a value not a percentile
%Process optional inputs
opts = ParseOptionalInputs(opts,varargin);

%load mask;
mask = load(opts.MaskDir);
mask = imresize(mask.brainoutline,size(frame));

%this is temporary additional mask to match the position of the 2022 data
temp = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\mask_posterioredge.mat');
mask = mask-temp.mask; mask(mask<0)=0;

%break mask into l and r hemisphere
temp = mask;
temp(:,34:end) = 0;
hemi_mask{1} = temp; 
temp = mask;
temp(:,1:34) = 0;
hemi_mask{2} = temp; 

nX = size(frame,1);
nY = size(frame,1);
    
%Bregma coordinates
bX = nX/opts.BregmaX+1;
bY = nY/opts.BregmaY;
    
% frame = SpatialGaussian(frame,3,1); %uncomment if you want to smooth
if opts.caxis_flag == 0
    climits = [prctile(frame(:),opts.caxis(1)),prctile(frame(:),opts.caxis(2))];
else
    climits = opts.caxis;
end
fH = cell(1,2);            
%loop through each hemisphere
for hemi = 1:2
    %Black background (to smooth the pixelated edges)
    cc = bwconncomp(hemi_mask{hemi},8);
    s = regionprops(cc,'Area','Centroid','ConvexHull');              
    fill(s(1).ConvexHull(:,1),s(1).ConvexHull(:,2),'k','LineWidth',2);

    cur_img = frame;
    cur_img(~hemi_mask{hemi})=0;


    %Plot the raw data  
    fH{hemi} = imagesc(cur_img); hold on

    %Make background transparent
    set(fH{hemi}, 'AlphaData',hemi_mask{hemi});            

    %Draw a smooth border to smooth the pixelations at the edges
    plot(s(1).ConvexHull(:,1),s(1).ConvexHull(:,2),'Color','w','LineWidth',3);

    %Set figure parameters
    colormap magma
    axis square
    caxis(climits)
    axis off
end %Hemi loop    
%Addbregma
scatter(bX,bY,600,'.','MarkerFaceColor',[1 0.2 .2],'MarkerEdgeColor',[1 0.2 .2]); 
ylim([0 68]);
xlim([0 68]);

%Set Axes
ylim([0 68]);
xlim([0 68]);
drawnow      

end %function end



