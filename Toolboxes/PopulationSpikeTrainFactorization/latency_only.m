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

function X1st = latency_only(X,len_trial)
% For each trial and each neuron, keep the first spike only
% Input arguments:
%  X         - Input matrix (size #cells x #time_points_per_trial * #trials)
%              composed of horizontally concatenated data trials
%  len_trial - Length of a trial in bins
% Output:
%  X1st      - Data matrix with first spikes only (same size as X)

[N,T] = size(X);
X1st = zeros(N,T);
n_trials = T/len_trial;
% Copy only 1st spike to X1st
for t = 1:n_trials
    for i = 1:N
        ind = find(X(i,(t-1)*len_trial+(1:len_trial)),1);
        if ~isempty(ind)
            X1st(i,(t-1)*len_trial+ind) = 1;
        end
    end
end

end
