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

function X_sh = shuffle_time(X,len_trial)
% Shuffle spikes times within trials to destroy timing information.
% Input arguments:
%  X         - Input matrix (size #cells x #time_points_per_trial * #trials)
%              composed of horizontally concatenated data trials
%  len_trial - Length of a trial in bins
% Output:
%  X_sh      - Shuffled data trials (same size as X)

[N,T] = size(X);
n_trials = T / len_trial;
X_sh = zeros(size(X));
% Shuffle bins for each trial:
for t = 1:n_trials
    % Draw random permutation
    p = randperm(len_trial);
    % Shuffle bins according to p
    X_sh(:,(t-1)*len_trial+(1:len_trial)) = X(:,(t-1)*len_trial+p);
end

end

