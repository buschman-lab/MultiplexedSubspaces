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

function [Acal,Wst] = stnmf(X,n_m,n_e,Wst)
% This function performs spatiotemporal NMF to factorize the data matrix X.
% Input arguments:
%  X    - Input matrix (size #cells x #time_points_per_sample * #samples)
%         composed of horizontally concatenated data samples
%  n_m  - Number of spatiotemporal modules
%  n_e  - Number of samples
% Optional input arguments:
%  Wst  - Spatiotemporal modules
%         (size #cells x #time_points_per_sample x #modules);
%         if given, modules are not optimized
% Output:
%  Acal - Activation coefficients (size #modules x #samples)
%  Wst  - Spatiotemporal modules
%         (size #cells x #time_points_per_sample x #modules)

optimize_Wst = nargin<4 || isempty(Wst); % If false, do not optimize modules

[N,T] = size(X);
n_t = T / n_e; % Number of bins per trial
X_nmf = zeros(N*n_t,n_e);
for i = 1:n_e
    X_nmf(:,i) = reshape(X(:,((i-1)*n_t+1):(i*n_t)),N*n_t,1);
end
if optimize_Wst
    [W,Acal] = nmf(X_nmf,n_m);
    Wst = zeros(N,n_t,n_m);
    for i = 1:n_m
        Wst(:,:,i) = reshape(W(:,i),N,n_t);
    end
else
    W = zeros(N*n_t,n_m);
    for i = 1:n_m
        W(:,i) = reshape(Wst(:,:,i),N*n_t,1);
    end
    [~,Acal] = nmf(X_nmf,n_m,W);
end

end

