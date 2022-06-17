function PlotMDSNetwork(data, varargin)

% PlotMDSNetwork - Makes a plot of data
%

% Parse inputs
opts.Colors = {'red', 'green', 'blue', 'black', 'magenta', 'cyan', '#DDAA33'};
opts.Radius = 0.03;
opts.NodeEdgeColor = 'none';
opts.FontSize = 12;
opts.FontWeight = 'bold';
opts.DrawEdges = 1;
opts.EdgeScaling = 4;
opts.EdgeR2Threshold = NaN;
opts.EdgePValThreshold = NaN;
opts.EdgeRThreshold = NaN;
opts.EdgeDataThreshold = NaN;
opts.PlotEdge = 'R2'; %can be 'none', 'R2', 'R', or 'pval' -- ignored if not doing correlation
opts.LabelEdge = 1;
opts.FillCircle = 1;
opts.SourceAreaName = 'TEST';
opts.AreaNames = {};
opts.EdgeText = 1;

opts.CalcCorrelation = 1;
opts.N_areas = 8;
opts.N_dims = 4;
opts.IgnoreNaNs = 1;

% Parse optional inputs
while length(varargin) >= 2,
    opts.(varargin{1}) = varargin{2};
    varargin = varargin(3:end);
end

if opts.CalcCorrelation,
    % Remove any NaNs
    if any(isnan(data(:))) && opts.IgnoreNaNs,
        for i = 1:size(data, 2),
            for j = 1:size(data, 2),
                good_ind = ~isnan(data(:, i)) & ~isnan(data(:, j));
                if sum(good_ind) == 0,
                    rho(i,j) = NaN; pval(i,j) = NaN;
                else
                    [rho(i,j), pval(i,j)] = corr(data(good_ind, i), data(good_ind, j));
                end
            end% data loop
        end %data loop
    else
        [rho, pval] = corr(data);
    end
    rho2 = rho.^2;

    %Make sure the diagonals are all 1 or MDS will complain (sometimes they
    %aren't because of machine precision errors)
    rho2([1:size(rho2, 1)] + size(rho2, 1)*([1:size(rho2, 1)]-1)) = 1;

    % Do MDS on the dissimilarities
    [Y, stress, disparities] = mdscale(1 - rho2, 2);
else
    % Do MDS on the dissimilarities
    [Y, stress, disparities] = mdscale(data, 2);
end



% Create plot
figure;

% Draw connections first
if ~strcmpi(opts.PlotEdge, 'none'),
    for i = 1:((opts.N_areas - 1)*opts.N_dims),
        for j = (i+1):((opts.N_areas - 1)*opts.N_dims),
            % Check whether there should be an edge between these points
            draw_edge = 0;
            if ~opts.CalcCorrelation,
                draw_edge = (data(i,j) >= opts.EdgeDataThreshold);
            elseif ~isnan(opts.EdgePValThreshold),
                draw_edge = (pval(i,j) <= opts.EdgePValThreshold);
            elseif ~isnan(opts.EdgeR2Threshold),
                draw_edge = (rho2(i,j) >= opts.EdgeR2Threshold);
            elseif ~isnan(opts.EdgeRThreshold),
                draw_edge = (abs(rho(i,j)) >= opts.EdgeRThreshold);
            end % draw edge?

            % Calculate edge width
            if ~opts.CalcCorrelation,
                edge_width = opts.EdgeScaling*data(i,j);
            elseif strcmpi(opts.PlotEdge, 'R2'),
                edge_width = opts.EdgeScaling*rho2(i,j);
            elseif strcmpi(opts.PlotEdge, 'R'),
                edge_width = opts.EdgeScaling*fisherZ(abs(rho(i,j)));
            elseif strcmpi(opts.PlotEdge, 'pval'),
                edge_width = -opts.EdgeScaling*log10(pval(i,j));
            end

            % Draw the edge
            if draw_edge,
                plot(Y([i j], 1), Y([i j], 2), 'k-', 'LineWidth', edge_width);
                hold on;
                if opts.LabelEdge,
                    mid_point = mean(Y([i j], 1:2), 1);
                    label_str = sprintf('%3.2f', edge_width/opts.EdgeScaling);
                    text(mid_point(1), mid_point(2), label_str, ...
                        'Color', '#000000', 'FontSize', opts.FontSize-2, 'FontWeight', opts.FontWeight, ...
                        'Rotation', atan(diff(Y([i j], 2))./diff(Y([i j], 1)))/pi*360, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                end
            end % draw edge
        end % dimension loop
    end %area loop
end %plot edges?

% Add nodes
for i = 1:(opts.N_areas - 1),
    for j = 1:opts.N_dims,
        % Create label
        label_str = '';
        if ~isempty(opts.AreaNames),
            label_str = opts.AreaNames{i};
        end
        label_str = strcat(label_str, num2str(j));

        if ~opts.FillCircle,
            plot(Y((i-1)*opts.N_dims + j, 1) + opts.Radius*sin([0:pi/100:2*pi]), ...
                Y((i-1)*opts.N_dims + j, 2) + opts.Radius*cos([0:pi/100:2*pi]), 'Color', opts.Colors{i});
            hold on;            
            text(Y((i-1)*opts.N_dims + j, 1), Y((i-1)*opts.N_dims + j, 2), label_str, 'Color', opts.Colors{i}, 'FontSize', opts.FontSize, 'FontWeight', opts.FontWeight, 'HorizontalAlignment', 'center');
            hold on;
        else
            fill(Y((i-1)*opts.N_dims + j, 1) + opts.Radius*sin([0:pi/100:2*pi]), ...
                Y((i-1)*opts.N_dims + j, 2) + opts.Radius*cos([0:pi/100:2*pi]), 'white', 'FaceColor', opts.Colors{i}, 'EdgeColor', opts.NodeEdgeColor);
            text(Y((i-1)*opts.N_dims + j, 1), Y((i-1)*opts.N_dims + j, 2), label_str, 'Color', 'white', 'FontSize', opts.FontSize, 'FontWeight', opts.FontWeight, 'HorizontalAlignment', 'center');
            hold on;
        end
    end
end
v = axis;
set(gca, 'XLim', [min(v(:)) max(v(:))], 'YLim', [min(v(:)) max(v(:))]);
axis square;
title(sprintf('Correlation of Projections from %s to Different Subspaces', opts.SourceAreaName));