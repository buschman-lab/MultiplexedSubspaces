function [hexagons] = MakeHexagonalGrid(surface_size, center_offset, hexagon_size)

% Makes a grid of hexagons that tiles a 2D surface. Returns as a list of
% polygons.
%
% Inputs:
%   surface_size - size of the 2D surface specified in 2 element array [width height]
%   center_offset - offset of the center hexagon. The tile goes in both
%                   directions from this starting point
%   hexagon_size - the length of each side of the regular hexagon