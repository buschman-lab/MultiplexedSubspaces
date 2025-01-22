function hex_h = PlotHexagonGridHistogram(ax_h, hex_count, pts, hexagons, hexagon_centers, varargin),

%Creates a plot of the hexagon histograms. 
%
% Inputs:
%  ax_h         -   Handle of the axis on which to plot. If empty or NaN,
%                   then a new axis is created. 
%  hex_count    -   A [M, D] matrix of the histogramed data for each
%                   of the M hexagons. D is the number of dimensions to
%                   plot. Each dimension is an independent pie wedge.
%  pts          -   A [Q, 2] matrix of the hexagon vertices.
%  hexagons     -   A [M, 6] matrix that contains the vertices of the M
%                   hexagons. The indices reference the rows of pts.

% Parse inputs
opts.NormalizeCountPerDimension = 0; % if true, then normalizes the counts to the maximum of a dimension across all hexagons
opts.NormalizeCountPerHexagon = 1; % if true, then normalizes counts within a hexagon
opts.NormalizeType = 'Max'; % can be 'Max', 'MinMax', or 'Sum' 

opts.PlotZ = NaN; % if a number, then plots using plot3 at that z height. Useful for plotting over images.
opts.HistogramEdgeColor = 'black'; % color of the edge
opts.HistogramEdgeStyle = 'none'; % transparency of the edge
opts.FaceColor = {[0.7 0.1 0.1], [0 0.25 0.5], [0.5 0.7 0.3], [0.6 0.2 0.5], [0.9 0.5 0.1], [0.6 0.5 0.55], [0.5 0.85 0.9], [0.2 0.35 0.38]};
opts.MaxFaceAlpha = 1;

opts.HexagonEdgeColor = 'black'; % color of the edge
opts.HexagonEdgeStyle = '-'; % transparency of the edge
opts.HexagonEdgeWidth = 2;

opts.Verbose = 1;

% Parse optional inputs
while length(varargin) >= 2,
    opts.(varargin{1}) = varargin{2};
    varargin = varargin(3:end);
end

if ~strcmpi(opts.NormalizeType, {'Max', 'Sum', 'MinMax'}),
    error('Unknown normalization type passed. Should be ''Max'', ''Sum'', or ''MinMax''');
end

% Check our axis
if isempty(ax_h) || ~ishandle(ax_h),
    fig_h = figure;
    ax_h = gca;
end

% Get the number of dimensions
M = size(hex_count, 1);
D = size(hex_count, 2);
if opts.Verbose,
    fprintf('Plotting hexagon grid with %d hexagons and %d dimensions.\n', M, D);
end

% Calculate the patches for inside one hexagon; we can just re-use the by
% recentering
cur_hex = 1;
cur_center = hexagon_centers(cur_hex, :); %center of current hexagon
D_patch = {};
for cur_D = 1:D,
    % Get intersection points of our angles
    [start_int_pt, start_int_section] = GetHexIntercept(pts(hexagons(cur_hex, :), :), cur_center, (cur_D - 1)/D*2*pi, 'Verbose', 0,'CurDim',cur_D);
    [end_int_pt, end_int_section] = GetHexIntercept(pts(hexagons(cur_hex, :), :), cur_center, cur_D/D*2*pi, 'Verbose', 0,'CurDim',cur_D);
    % now store the pts for the section patch
    D_patch{cur_D} = cat(1, cur_center, start_int_pt, pts(hexagons(cur_hex, (mod(start_int_section, 6)+1):end_int_section), :), end_int_pt);
    % Make everything relative to center for easier addition...
    D_patch{cur_D} = D_patch{cur_D} - repmat(cur_center, [size(D_patch{cur_D}, 1) 1]);
end % dimension loop
if opts.Verbose,
    fprintf('Created patches.\n');
end

% Normalize counts across hexagons?
if opts.NormalizeCountPerDimension,
    if strcmpi(opts.NormalizeType, 'Sum'),
        hex_count = hex_count./repmat(sum(hex_count, 1), [size(hex_count, 1) 1]);
    elseif any(strcmpi(opts.NormalizeType, {'Max', 'MinMax'})),
        if strcmpi(opts.NormalizeType, 'MinMax'),
            hex_count = hex_count - repmat(min(hex_count, [], 1), [size(hex_count, 1) 1]);
        end
        hex_count = hex_count./repmat(max(hex_count, [], 1), [size(hex_count, 1) 1]);
    end
    if opts.Verbose,
        fprintf('Normalized counts by dimension.\n');
    end
end

% Normalize the counts within a hexagon
if opts.NormalizeCountPerHexagon,
    if strcmpi(opts.NormalizeType, 'Sum'),
        hex_count = hex_count./repmat(sum(hex_count, 2), [1 size(hex_count, 2)]);
    elseif any(strcmpi(opts.NormalizeType, {'Max', 'MinMax'})),
        if strcmpi(opts.NormalizeType, 'MinMax'),
            hex_count = hex_count - repmat(min(hex_count, [], 2), [1 size(hex_count, 2)]);
        end
        hex_count = hex_count./repmat(max(hex_count, [], 2), [1 size(hex_count, 2)]);
    end        
    if opts.Verbose,
        fprintf('Normalized counts per hexagon.\n');
    end
end


% Loop through the hexagons, plotting as we go
hex_h = [];
axes(ax_h);
if opts.Verbose,
    fprintf('Creating hexagons...\n');
end
for cur_hex = 1:M,
    % Create the patches for each dimension
    for cur_D = 1:D,
        cur_color = opts.FaceColor{mod(cur_D - 1, length(opts.FaceColor)) + 1};
        cur_alpha = opts.MaxFaceAlpha*hex_count(cur_hex, cur_D);
        if ~isnan(opts.PlotZ),
            cur_h = patch('XData', D_patch{cur_D}(:, 1) + hexagon_centers(cur_hex, 1), ...
                'YData', D_patch{cur_D}(:, 2) + hexagon_centers(cur_hex, 2), ...
                'ZData', opts.PlotZ*ones(size(D_patch{cur_D}, 1), 1), ...
                'FaceColor', cur_color, 'FaceAlpha', cur_alpha, ...
                'EdgeColor', opts.HistogramEdgeColor, 'LineStyle', opts.HistogramEdgeStyle);
        else
            cur_h = patch('XData', D_patch{cur_D}(:, 1) + hexagon_centers(cur_hex, 1), ...
                'YData', D_patch{cur_D}(:, 2) + hexagon_centers(cur_hex, 2), ...
                'FaceColor', cur_color, 'FaceAlpha', cur_alpha, ...
                'EdgeColor', opts.HistogramEdgeColor, 'LineStyle', opts.HistogramEdgeStyle);
        end
        hex_h = cat(1, hex_h, cur_h);
    end %dimension loop
    if opts.Verbose && (mod(cur_hex, floor(M/10)) == 0),
        fprintf('\tPlotted %3.1f%% of the hexagons...\n', cur_hex/M*100);
    end
end % hexagon loop
if opts.Verbose,
    fprintf('\tDone.\n');
end

% Plot the hexagons on top
if ~strcmpi(opts.HexagonEdgeStyle, 'none'),
    if opts.Verbose,
        fprintf('Adding hexagon grid on top.\n');
    end
    hex_h = cat(1, PlotHexagonGrid(ax_h, pts, hexagons, 'EdgeColor', opts.HexagonEdgeColor, 'EdgeStyle', opts.HexagonEdgeStyle, 'EdgeWidth', opts.HexagonEdgeWidth, 'FaceColor', 'none', 'FaceAlpha', 0));
end

end %PlotHexagonGridHistogram function



