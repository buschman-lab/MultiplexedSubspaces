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

function counts = count_spikes(trains,bin_size,len_trial)
% Bins the spike trains by counting the number of spikes in each bin.
% Input arguments:
%  trains    - Input cell array (size #cells x #trials)
%              with vectors of spike times in milliseconds
%  bin_size  - Length of a bin in milliseconds
%  len_trial - Length of a trial in milliseconds
% Output:
%  counts    - Cell array (length #trials) with spike count matrices of
%              size (#cells x #bins_per_trial)

n_bins = ceil(len_trial/bin_size);
[N,n_trials] = size(trains);
counts = cell(n_trials,1);
for l = 1:n_trials
    counts{l} = zeros(N,n_bins);
    for i = 1:N
        if n_bins > 0
            for j = 1:n_bins
                counts{l}(i,j) = sum((trains{i,l}>(j-1)*bin_size) & (trains{i,l} <= j*bin_size));
            end
        end
    end
end

end
