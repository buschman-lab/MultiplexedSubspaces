function Figure6b(rho,sig_thresh)
opts.Threshold = 0.5; %threshold for number of pixels within a hexagon in order to fill hex with a color. 
N_dims = 2:7; %dimensions to use for plotting (def 2:7)
fp = fig_params_cortdynamics;

%%
%threshold and binarize based on permutation-derived significance threshold
for i = 1:10
   temp = rho(:,:,i);
   temp(abs(temp)<sig_thresh(i))=0;
   rho(:,:,i) = temp;
end
%Binarize. 
rho(rho>0)=1; 
rho(rho<=0)=0;

%% Parse the hexagons
hex_data = rho(:,:,N_dims);
% hex_data = rho(:,:,5);

data_x = linspace(0, 68, 68);
data_y = linspace(0, 68, 68);

% Create our hexagon grid over an imaginary 100x100 pixel space
%position and radius are the best parameteres to adjust here
%good sizes include 5 with center pos 5,4.5 and 6 with centerpos 7,0
[pts, hexagons, hexagon_centers] = CreateHexagonGrid('RemoveDuplicateVertices', 1, 'Radius', 5, 'GridSize', [0 100; 100 0], 'CenterPos', [5 4.5], 'RotationAngle', 0, 'Verbose', 0);

% Histogram it
hex_count = HexagonGridHistogram(hex_data, data_x, data_y, pts, hexagons);
% Remove empty hexagons
badhex = sum(isnan(hex_count),2)==size(hex_count,2);
hex_count(badhex,:)=[];
hexagon_centers(badhex,:)=[];
hexagons(badhex,:)=[];

%Replaces empty sections with zeros
hex_count(isnan(hex_count))=0;

% Plot it
hex_count(hex_count<opts.Threshold)=0;
hex_count(hex_count>=opts.Threshold)=1;
figure;
hex_h = PlotHexagonGridHistogram(gca, hex_count, pts, hexagons, hexagon_centers, 'NormalizeCountPerDimension', 0,... %no normalization needed (will through errs since division by zero
    'NormalizeCountPerHexagon', 0, 'Verbose', 0,'HexagonEdgeColor',[0.25 0.25 0.25],'HexagonEdgeWidth',3.5);
set(gca,'ydir','reverse')
axis off            

title('Frontal Motor','fontweight','normal');
fp.FigureSizing(gcf,[1 1 8 8],[10 10 10 10])
end

