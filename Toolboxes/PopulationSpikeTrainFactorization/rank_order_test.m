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

function groups = rank_order_test(X,len_trial,templates)
% Estimates predictions for each group based on rank order templates,
% following Panzeri and Diamond, 2010.
% Input arguments:
%  X         - Input matrix (size #cells x #time_points_per_trial * #trials)
%              composed of horizontally concatenated data samples
%  len_trial - Length of a trial in bins
%  templates - Templates obtained from rank_order_train
% Output:
%  groups    - Class label predictions for each trial
%              (vector of length #trials)

[N,T] = size(X);
n_trials = T/len_trial;
groups = zeros(n_trials,1);
% Compute first spike ranks of all trials
for t = 1:n_trials
    latencies = ones(N,1)*(len_trial+1);
    for i = 1:N
        ind = find(X(i,(t-1)*len_trial+(1:len_trial)),1);
        if ~isempty(ind)
            latencies(i) = ind;
        end
    end
    rho = corr(tiedrank(latencies),templates,'type','Spearman');
    [~,groups(t)] = max(rho);
end

end

