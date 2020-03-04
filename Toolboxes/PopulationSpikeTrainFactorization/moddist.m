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

function [d,order,signs] = moddist(m1,m2,...
                                   noneg,geodesic)
% Calculate module distance based on geodesic distance for best matching
% modules. Modules are assumed to be column vectors. 
% Input arguments:
%  m1       - First matrix of column modules
%  m2       - Second matrix of column modules
% Optional input arguments:
%  noneg    - By default module and minus module are considered to have
%             zero distance; Set noneg to true if they should be
%             considered different
%  geodesic - Use geodesic distance if true and euclidean distance if
%             false
% Output:
%  d        - Distance between module sets
%  order    - Closest permutation found of second module set
%  signs    - Closest signs of second module set

if nargin<3 || isempty(noneg)
    noneg = false;
end
if nargin<4 || isempty(geodesic)
    geodesic = true;
end
% Normalize vectors to unit length
m1 = m1 ./ repmat(sqrt(sum(m1.^2)),size(m1,1),1);
m2 = m2 ./ repmat(sqrt(sum(m2.^2)),size(m2,1),1);
n_m = size(m1,2);
% Greedily order modules to calculate minimal distance
dm = zeros(1,n_m);
% Also try negative
dmneg = zeros(1,n_m);
% Order of second module set
order = 1:n_m;
signs = ones(1,n_m);
for i = 1:n_m
    for j = i:n_m
        if geodesic
            % Geodesic distance
            dm(j) = acos(dot(m1(:,i),m2(:,order(j))))./pi;
            % Geodesic distance to negative modules
            dmneg(j) = acos(dot(m1(:,i),-m2(:,order(j))))./pi;
        else
            % Euclidean distance
            dm(j) = norm(m1(:,i)-m2(:,order(j)));
            % Euclidean distance to negative modules
            dmneg(j) = norm(m1(:,i) + m2(:,order(j)));
        end
    end
    if noneg || min(dm(i:n_m))<=min(dmneg(i:n_m))
        [dm(i),j_min] = min(dm(i:n_m));
    else
        [dm(i),j_min] = min(dmneg(i:n_m));
        signs(i) = -1;
    end
    t = order(j_min);
    order(j_min) = order(i);
    order(i) = t;
end
d = mean(dm);

end

