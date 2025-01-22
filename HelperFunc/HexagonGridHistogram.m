function hex_count = HexagonGridHistogram(data, data_x, data_y, pts, hexagons, varargin),

%Creates a histogram of data with hexagon grid. 
%
% Inputs:
%  data         -   A matrix of to-be-counted data. Should be [x, y, N] where
%                   x and y are the length of the data_x and data_y vectors
%                   and N is the number of dimensions to count.
%  data_x       -   A vector of the x-positions of the data matrix.
%  data_y       -   A vector of the y-positions of the data matrix.
%  pts          -   A [Q, 2] matrix of the hexagon vertices.
%  hexagons     -   A [M, 6] matrix that contains the vertices of the M
%                   hexagons. The indices reference the rows of pts.

% Parse inputs
opts.IgnoreNaN = 1; %ignore any NaNs

% Parse optional inputs
while length(varargin) >= 2,
    opts.(varargin{1}) = varargin{2};
    varargin = varargin(3:end);
end

% Check our inputs
if (length(data_x) ~= size(data, 1)) || (length(data_y) ~= size(data, 2)),
    error('The data x and y vectors must match the size of the data matrix.');
end

% Get the number of dimensions
N = size(data, 3);

% Tile the full space with the data_x and data_y vectors
[X, Y] = meshgrid(data_x(:), data_y(:));

% Loop through the hexagons, counting as we go
hex_count = zeros(size(hexagons, 1), N);
for cur_hex = 1:size(hexagons, 1),
    in = inpolygon(X(:), Y(:), pts(hexagons(cur_hex, :), 1), pts(hexagons(cur_hex, :), 2));
    if any(in), 
        [in_row, in_col] = ind2sub(size(X), find(in));
        if opts.IgnoreNaN,
            hex_count(cur_hex, :) = squeeze(mean(data(in_row, in_col, :), [1 2], 'omitnan'));
        else
            hex_count(cur_hex, :) = squeeze(mean(data(in_row, in_col, :), [1 2], 'includenan'));
        end
    else
        hex_count(cur_hex, :) = NaN;
    end
end % hexagon loop
