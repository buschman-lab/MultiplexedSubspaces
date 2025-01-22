function CorticalHexPlots(data_all,cur_rec,motif,area,varargin)
%Camden - timeless
%% parse optional inputs
opts.Threshold = 0.5; %percent of pixels to be considered significant, only if opts.Binarize==1
opts.Binarize = 1; %if 1 (def) then plot the whether a region is y/n involved. If 0 then use the actual correlation value
opts.N_dims = 2:7; %dimensions to use for plotting (def 2:7)
opts = ParseOptionalInputs(opts,varargin); 

fp = fig_params_cortdynamics;
%% Load the data
data = data_all{cur_rec}; %grab single recording
%get the motif and area index
idx = find(cat(1,data.cur_motif)==motif & cat(1,data.cur_a)==area);

rho_all = data(idx).rho_all; 
%direction of correlation is arbitrary. Flip so strongest are positive -
%that is our subspace network for this dimension. 
for i = 1:size(rho_all,3)
    if nansum(rho_all(:,:,i),'all') < nansum(-1*rho_all(:,:,i),'all')
        rho_all(:,:,i)=-1*rho_all(:,:,i);
    end
end

if opts.Binarize ==1
    %parse data structure
    rho = rho_all; 
    sig_thresh = data(idx).sig_thresh; 
    %threshold and binarize
    for i = 1:10
       temp = rho(:,:,i);
       temp(abs(temp)<sig_thresh(i))=0;
       rho(:,:,i) = temp;
    end
    %binarize. anything that is significantly positively correlated is our network
    rho(rho>0)=1; 
    rho(rho<=0)=0;

elseif opts.Binarize==0
    rho = rho_all; %... there can be weak negative correlations here (so anticorrelated with the subspace?) ... exclude those right now, because we want to see the positively correlated network
    rho(rho<0)=0; 
else
    error('unknown Binarization option (must be zero or one');
end

% to get the maximum correlation with other areas, uncomment
% a = rho(:,:,6);
% b = rho(:,:,9);
% c = a(:)+b(:);
% ovr=(2*sum(c==2))/(sum(a(:)==1)+sum(b(:)==1));
% spatialcorr=corr(a(:),b(:),'rows','complete');
% %maximum correlation
% a = rho(:,:,10);
% r = NaN(1,10);
% rr =NaN(1,10);
% for i = 1:10
%     b = rho(:,:,i);
%     r(i) = corr(a(:),b(:),'rows','complete');
%     c = a(:)+b(:);
%     rr(i) = (2*sum(c==2))/(sum(a(:)==1)+sum(b(:)==1));
% end
% 
% %


%% Parse the hexagons
hex_data = rho(:,:,opts.N_dims);
% hex_data = rho(:,:,5);

data_x = linspace(0, 68, 68);
data_y = linspace(0, 68, 68);


% Create our hexagon grid over an imaginary 100x100 pixel space
%position and radius are the best parameteres to adjust here
%good sizes include 5 with center pos 5,4.5 and 6 with centerpos 7,0
% [pts, hexagons, hexagon_centers] = CreateHexagonGrid('RemoveDuplicateVertices', 1, 'Radius', 6, 'GridSize', [0 100; 100 0], 'CenterPos', [7 0], 'RotationAngle', 0, 'Verbose', 0);
[pts, hexagons, hexagon_centers] = CreateHexagonGrid('RemoveDuplicateVertices', 1, 'Radius', 5, 'GridSize', [0 100; 100 0], 'CenterPos', [5 4.5], 'RotationAngle', 0, 'Verbose', 0);

% Histogram it
hex_count = HexagonGridHistogram(hex_data, data_x, data_y, pts, hexagons);
% Remove empty hexagons
if opts.Binarize==1
    badhex = sum(isnan(hex_count),2)==size(hex_count,2); %just remove those outside the brain
else
    badhex = sum(isnan(hex_count),2)==size(hex_count,2) | sum(hex_count==0,2)==size(hex_count,2); %also remove empty to avoid normalization error
end
hex_count(badhex,:)=[];
hexagon_centers(badhex,:)=[];
hexagons(badhex,:)=[];

%Replaces empty sections with zeros
hex_count(isnan(hex_count))=0;

% Plot it
if opts.Binarize ==1
    hex_count(hex_count<opts.Threshold)=0;
    hex_count(hex_count>=opts.Threshold)=1;
    figure;
    hex_h = PlotHexagonGridHistogram(gca, hex_count, pts, hexagons, hexagon_centers, 'NormalizeCountPerDimension', 0,... %no normalization needed (will through errs since division by zero
        'NormalizeCountPerHexagon', 0, 'Verbose', 0,'HexagonEdgeColor',[0.25 0.25 0.25],'HexagonEdgeWidth',3.5);
    set(gca,'ydir','reverse')
    axis off            
else
    figure;
    hex_h = PlotHexagonGridHistogram(gca, hex_count, pts, hexagons, hexagon_centers, 'NormalizeCountPerDimension', 0,...
        'NormalizeCountPerHexagon', 1, 'NormalizeType', 'MinMax', 'Verbose', 0,'HexagonEdgeColor',[0.25 0.25 0.25],'HexagonEdgeWidth',3.5);
    set(gca,'ydir','reverse')
    axis off        
end

end %function end














