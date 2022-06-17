close all; clear all;

% Create our hexagon grid over an imaginary 100x100 pixel space
[pts, hexagons, hexagon_centers] = CreateHexagonGrid('RemoveDuplicateVertices', 1, 'Radius', 10, 'GridSize', [-100 100; 100 -100], 'CenterPos', [0 0], 'RotationAngle', 0, 'Verbose', 0);

% We can plot the hexagon grid to make sure it looks okay
figure; PlotHexagonGrid(gca, pts, hexagons); axis equal;

%% Now create some fake data with sparse noise
N_dims = 4;
data_x = linspace(-100, 100, 201);
data_y = linspace(-100, 100, 201);
data = rand([length(data_x), length(data_y), N_dims]) >= 0.9;
% Smooth it
sm_kern = ones(5, 5, 1); sm_kern = sm_kern./sum(sm_kern, "all");
data = convn(data, sm_kern, 'same');

% Visualize the random data
figure;
for i = 1:N_dims
    subplot(floor(sqrt(N_dims)), ceil(sqrt(N_dims)), i);
    subtitle(sprintf('Dimension %d', i));
    imagesc(data(:, :, i));
end

%% Histogram it
hex_count = HexagonGridHistogram(data, data_x, data_y, pts, hexagons);

%% Now plot it!
%figure;
%hex_h = PlotHexagonGridHistogram(gca, hex_count, pts, hexagons, hexagon_centers, 'NormalizeCountPerDimension', 0, 'NormalizeCountPerHexagon', 1, 'NormalizeType', 'Max', 'Verbose', 1);

figure;
hex_h = PlotHexagonGridHistogram(gca, hex_count, pts, hexagons, hexagon_centers, 'NormalizeCountPerDimension', 0, 'NormalizeCountPerHexagon', 1, 'NormalizeType', 'MinMax', 'HexagonEdgeWidth', 2, 'Verbose', 1);