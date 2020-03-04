% Copyright (C) 2016 Arno Onken
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <http://www.gnu.org/licenses/>.

function o = overlap(X)
% Calculates overlap as the geodesic similarity averaged over all columns.
% Input argument:
%  X - Input data matrix; each column of X is assumed to be a vector
% Output:
%  o - Mean overlap (between 0  for no overlap and 1 for maximum overlap)

% Normalize vectors to unit length
X = X ./ repmat(sqrt(sum(X.^2)),size(X,1),1);
% Average overlap over all pairs
o = 0;
for i = 2:size(X,2)
    for j = 1:(i-1)
        % Geodesic similarity
        o = o + (1 - 2 * acos(dot(X(:,i),X(:,j)))./pi);
    end
end
o = o ./ (size(X,2)*(size(X,2)-1)/2);

end

