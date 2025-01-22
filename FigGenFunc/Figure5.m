function Figure5(clust_colmap,clust_size,k,sorted_sim_mat)

% Process optional inputs
opts.ColorBarLabel_Offset = 0;
opts.ColorBarLabel_Height = 0.025;

% Create figure with the similarity between clusters
clust_sim_fig = figure;
imagesc(sorted_sim_mat);
hold on;
axis square;
axis off;
v = axis;
title('Similarity of Subspace Networks','Fontweight','normal')

% Create rectangles on the outside to indicate the clusters
rect_edges = [0; cumsum(clust_size)];
for i = 1:k,
    % Along bottom
    h = rectangle('Position', [(rect_edges(i)+0.5) (1+opts.ColorBarLabel_Offset)*rect_edges(end) clust_size(i) opts.ColorBarLabel_Height*rect_edges(end)], 'FaceColor', clust_colmap(i, :), 'EdgeColor', 'none');
    set(h, 'Clipping', 'off');
    text((rect_edges(i)+clust_size(i)/2+0.5), (1+opts.ColorBarLabel_Offset + opts.ColorBarLabel_Height/2)*rect_edges(end), sprintf('%d', i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'white');
    % Along left side
    h = rectangle('Position', [-(opts.ColorBarLabel_Offset+opts.ColorBarLabel_Height)*rect_edges(end) (rect_edges(i)+0.5) opts.ColorBarLabel_Height*rect_edges(end) clust_size(i)], 'FaceColor', clust_colmap(i, :), 'EdgeColor', 'none');
    set(h, 'Clipping', 'off');
    text(-(opts.ColorBarLabel_Offset+opts.ColorBarLabel_Height/2)*rect_edges(end), (rect_edges(i)+clust_size(i)/2+0.5), sprintf('%d', i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'white');
    % Draw white box
    plot([rect_edges(i) rect_edges(i) rect_edges(i+1) rect_edges(i+1) rect_edges(i)], [rect_edges(i) rect_edges(i+1) rect_edges(i+1) rect_edges(i) rect_edges(i)], 'w-', 'LineWidth', 2);
end


end

