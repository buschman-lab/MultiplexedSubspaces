function [pts, hexagons, hexagon_centers] = CreateHexagonGrid(varargin),

%Creates a grid of hexagons in 2D. Returns a list of points that are the
%vertices of the hexagons as well as an array of grouping of the index

% Parse inputs
opts.CenterPos = [0 0]; % center position of the hexagon grid
opts.Radius = 10; % radius (in points) of the hexagon
opts.RotationAngle = 0; % rotation (in degrees) of the grid
opts.GridSize = [-100 100; 100 -100]; % the shape of the rectangle that 
                                     % should be tiled. Only hexagons with 
                                     % ALL points within the grid will be kept.
                                     % Grid size should be in the format of 
                                     % upper-left corner in first row and 
                                     % lower-right in second row 
opts.RemoveDuplicateVertices = 0;
opts.UniqueThreshold = 10^-5;
opts.TrimPartial = 1; %whether to trim hexagons that partially extend beyond borders
opts.Verbose = 0;

% Parse optional inputs
while length(varargin) >= 2,
    opts.(varargin{1}) = varargin{2};
    varargin = varargin(3:end);
end

if opts.Verbose, fprintf('Creating hexagon grid with hexagon radius (side length) of %3.1f, centered on [%3.1f, %3.1f], rotated %3.1f degrees.\n', opts.Radius, opts.CenterPos, opts.RotationAngle); end

% Create a function for rotation matrices
rotMat = @(ang) ([cos(ang) -sin(ang); sin(ang) cos(ang)]);
% Convert input to radians
opts.RotationAngle = opts.RotationAngle/180*pi;

% To make sure we can rotate the grid, expand it out
cur_grid_size = opts.GridSize;
cur_grid_size(1, 1) = cur_grid_size(1, 1) - (opts.GridSize(2, 1) - opts.GridSize(1, 1))/2;
cur_grid_size(2, 1) = cur_grid_size(2, 1) + (opts.GridSize(2, 1) - opts.GridSize(1, 1))/2;
cur_grid_size(1, 2) = cur_grid_size(1, 2) + (opts.GridSize(1, 2) - opts.GridSize(2, 2))/2;
cur_grid_size(2, 2) = cur_grid_size(2, 2) - (opts.GridSize(1, 2) - opts.GridSize(2, 2))/2;

% Make sure that rectangle is valid
if (opts.GridSize(2,1) <= opts.GridSize(1,1)) || (opts.GridSize(2,2) >= opts.GridSize(1,2)),
    error('GridSize must be specified from upper-left corner to lower-right corner.');
end
if opts.Verbose,
    fprintf('Limiting grid to rectangle from [%3.1f, %3.1f] to [%3.1f, %3.1f].\n', opts.GridSize);
end

% Start with a vector that will form the x-center of the hexagon grid
if (opts.CenterPos(1) <= cur_grid_size(2, 1)) & (opts.CenterPos(1) >= cur_grid_size(1,1)),
    % Center is within the grid, so go in both directions
    base_hexagon_centers = [opts.CenterPos(1):(3*opts.Radius):cur_grid_size(2, 1)];
    base_hexagon_centers = cat(2, fliplr([opts.CenterPos(1):(-3*opts.Radius):cur_grid_size(1,1)]), base_hexagon_centers(2:end));
elseif (opts.CenterPos(1) <= cur_grid_size(1, 1)),
    % Center is to the left of the grid, so only go to the right
    base_hexagon_centers = [opts.CenterPos(1):(3*opts.Radius):cur_grid_size(2, 1)];
    base_hexagon_centers = base_hexagon_centers(base_hexagon_centers >= cur_grid_size(1,1));
elseif (opts.CenterPos(1) >= cur_grid_size(2, 1)),
    % Center is to the right of the grid, so only go to the left
    base_hexagon_centers = fliplr([opts.CenterPos(1):(-3*opts.Radius):cur_grid_size(1,1)]);
    base_hexagon_centers = base_hexagon_centers(base_hexagon_centers <= cur_grid_size(2,1));
end
% Add second dimension
base_hexagon_centers = cat(1, base_hexagon_centers, opts.CenterPos(2)*ones(1, size(base_hexagon_centers, 2)));
base_hexagon_centers = base_hexagon_centers';
if opts.Verbose,
    fprintf('Baseline strip of hexagons has points: %s\n', mat2str(base_hexagon_centers));
end

% Do the same, but now for the y-center of the hexagon grid
y_step = sin(pi/3)*opts.Radius; 
x_step = 1.5*opts.Radius;
cur_y = 0; 
cur_x = 0;
hexagon_centers = [];
if (opts.CenterPos(2) <= cur_grid_size(1, 2)) && (opts.CenterPos(2) >= cur_grid_size(2,2)),
    % Center is within the grid, so go in both directions
    
    % Go up to the top of the grid
    while (cur_y <= cur_grid_size(1, 2)),
        hexagon_centers = cat(1, hexagon_centers, base_hexagon_centers + [cur_x cur_y]);
        cur_y = cur_y + y_step;
        cur_x = cur_x + x_step;
        x_step = -x_step;
    end
    % Go down to the bottom of the grid
    x_step = 1.5*opts.Radius;
    cur_y = 0;
    cur_x = 0;
    while (cur_y >= cur_grid_size(2, 2)),
        hexagon_centers = cat(1, hexagon_centers, base_hexagon_centers + [cur_x cur_y]);
        cur_y = cur_y - y_step;
        cur_x = cur_x + x_step;
        x_step = -x_step;
    end
elseif (opts.CenterPos(2) <= cur_grid_size(2, 2)),
    % Center is below the grid, so only go up

    % Go up to the top of the grid
    while (cur_y <= cur_grid_size(1, 2)),
        hexagon_centers = cat(1, hexagon_centers, base_hexagon_centers + [cur_x cur_y]);
        cur_y = cur_y + y_step;
        cur_x = cur_x + x_step;
        x_step = -x_step;
    end

elseif (opts.CenterPos(1) >= cur_grid_size(2, 1)),
    % Center is above the grid, so only go down 
    
    % Go down to the bottom of the grid
    cur_y = 0;
    cur_x = 0;
    while (cur_y >= cur_grid_size(2, 2)),
        hexagon_centers = cat(1, hexagon_centers, base_hexagon_centers + [cur_x cur_y]);
        cur_y = cur_y - y_step;
        cur_x = cur_x + x_step;
        x_step = -x_step;
    end
end

% Now generate all of the vertices of the hexagon grid
if opts.Verbose, fprintf('Generating all of the hexagon vertices...\n'); end
base_hex_pts = [cos([0:(pi/3):(5*pi/3)]') sin([0:(pi/3):(5*pi/3)]')]*opts.Radius;
pts = [];
hexagons = NaN*ones(size(hexagon_centers, 1), 6);
% Loop through all of the centers, copying
for cur_center = 1:size(hexagon_centers, 1),
    pts = cat(1, pts, hexagon_centers(cur_center, :) + base_hex_pts);
    hexagons(cur_center, :) = (cur_center-1)*6 + [1:6];
end

% Remove any duplicates?
if opts.RemoveDuplicateVertices,
    if opts.Verbose, fprintf('Removing duplicate vertices...\n'); end
    uniq_pts = pts(1, :);
    for cur_center = 1:size(hexagon_centers, 1),
        % Loop through all of the points
        for i = 1:size(hexagons, 2),
            % Is this in our unique list?
            uniq_pts_dist = sqrt(sum((repmat(pts(hexagons(cur_center, i), :), [size(uniq_pts, 1) 1]) - uniq_pts).^2, 2));
            ind = find(uniq_pts_dist <= opts.UniqueThreshold, 1, 'first');
            if ~isempty(ind),
                % Already in unique list
                hexagons(cur_center, i) = ind;
            else
                % Not in unique list
                uniq_pts = cat(1, uniq_pts, pts(hexagons(cur_center, i), :));
                hexagons(cur_center, i) = size(uniq_pts, 1);
            end
        end % hexagon point loop
    end % hexagon loop
    pts = uniq_pts;
end

if opts.Verbose, fprintf('Transforming hexagons...\n'); end

% Rotate all of the points
pts = (rotMat(opts.RotationAngle)*pts')';

% Trim hexagons
if opts.Verbose, fprintf('Trimming hexagons outside the range of the grid...\n'); end
% Create a list of bad points (outside the grid size)
bad_ind = (pts(:, 1) < opts.GridSize(1,1)) | (pts(:, 1) > opts.GridSize(2,1)) | (pts(:, 2) > opts.GridSize(1,2)) | (pts(:, 2) < opts.GridSize(2,2));
if opts.Verbose,
    figure; plot(pts(:, 1), pts(:, 2), 'bo');
    hold all;
    plot(pts(bad_ind, 1), pts(bad_ind, 2), 'rx');
end
if opts.TrimPartial,
    % Trim anything that is outside the window    
    bad_hex = any(ismember(hexagons, find(bad_ind)), 2);
else
    % Only trim those hexagons fully outside the window
    bad_hex = all(ismember(hexagons, find(bad_ind)), 2);
end
% Take out any bad hexagons
hexagons = hexagons(~bad_hex, :);
hexagon_centers = hexagon_centers(~bad_hex, :);
% Take out any points that aren't reference by hexagons anywhere
used_pts = ismember([1:size(pts, 1)], unique(hexagons(:)));
pt_ind = [1:size(pts, 1)];
pt_ind(used_pts) = [1:sum(used_pts)];
hexagons = pt_ind(hexagons);
pts = pts(used_pts, :);

if opts.Verbose, 
    fprintf('Finished creating hexagon grid with %d hexagons.\n', size(hexagons, 1));
    figure;
    PlotHexagonGrid(gca, pts, hexagons); axis equal;
end

end
