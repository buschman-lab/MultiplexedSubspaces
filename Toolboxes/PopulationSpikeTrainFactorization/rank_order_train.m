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

function templates = rank_order_train(X,len_trial,groups)
% Calculates mean rank templates for each group, following Panzeri and
% Diamond, 2010.
% Input arguments:
%  X         - Input matrix (size #cells x #time_points_per_trial * #trials)
%              composed of horizontally concatenated data samples
%  len_trial - Length of a trial in bins
%  groups    - Class labels of each trial (e.g. stimulus identity);
%              (vector of length #trials)
% Output:
%  templates - Class templates for classifying

[N,T] = size(X);
n_trials = T / len_trial;
id = min(groups):max(groups);
ranks = cell(length(id),1);
for i = 1:length(id)
    ranks{i} = zeros(N,sum(groups==id(i)));
end
g_ind = ones(length(id),1);
% Compute first spike ranks of all trials
for t = 1:n_trials
    latencies = ones(N,1)*(len_trial + 1);
    for i = 1:N
        ind = find(X(i,(t-1)*len_trial+(1:len_trial)),1);
        if ~isempty(ind)
            latencies(i) = ind;
        end
    end
    ranks{groups(t)}(:,g_ind(groups(t))) = tiedrank(latencies);
    g_ind(groups(t)) = g_ind(groups(t)) + 1;
end
% Calculate average ranks
templates = zeros(N,length(id));
for i = 1:length(id)
    templates(:,i) = tiedrank(mean(ranks{i},2));
end

end

