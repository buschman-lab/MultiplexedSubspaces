% Get point where two lines intercept
function [int_pt, int_section] = GetHexIntercept(hex_pts, origin, cur_ang, varargin),

opts.Verbose = 1; %should be verbose?
opts.CurDim = 0; %include for plotting if verbose

% Parse optional inputs
while length(varargin) >= 2,
    opts.(varargin{1}) = varargin{2};
    varargin = varargin(3:end);
end



% Get current points
cur_pts = hex_pts - origin;

%Rotate by our angle
rotMat = @(ang) ([cos(ang) -sin(ang); sin(ang) cos(ang)]);
% Convert input to radians
cur_pts = (rotMat(-cur_ang)*cur_pts')';

int_pt = NaN*ones(1, 2);
for cur_section = 1:6,
    section_pts = cur_pts([cur_section mod(cur_section, 6) + 1], :);
    % Does this span the x axis?
    if (section_pts(1, 2) < 0) && (section_pts(2, 2) >= 0),
        % This is our section
        int_section = cur_section;
        % Find intersection
        int_pt(1) = section_pts(1, 1) - section_pts(1, 2)*diff(section_pts(:, 1))/diff(section_pts(:, 2));
        int_pt(2) = 0;
        %break; % no more searching
    end
end % loop through sections

%Unrotate
int_pt = (rotMat(cur_ang)*int_pt')' + origin;

if opts.Verbose,
    figure;
    plot(hex_pts(:, 1), hex_pts(:, 2), 'ro');
    hold all;
    plot(origin(1), origin(2), 'kx');
    plot(int_pt(1), int_pt(2), 'gd');
    len = 3;    
    title(sprintf('Cur Dim %d',opts.CurDim),'fontweight','normal'); 
    plot(origin(1) + [0 len*cos(cur_ang)], origin(2) + [0 len*sin(cur_ang)], 'g-');
end