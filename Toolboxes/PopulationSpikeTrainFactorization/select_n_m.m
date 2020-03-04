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

function n_m = select_n_m(X_train,groups_train,n_e_train,X_test,groups_test,n_e_test)
% Selects the optimal number of spatiotemporal modules.
% Input arguments:
%  X_train      - Training input matrix (size #cells
%                 x #time_points_per_trial * #training_set_trials)
%                 composed of horizontally concatenated population responses
%  groups_train - Class labels of each trial (e.g. stimulus identity);
%                 (vector of length #training_set_trials)
%  n_e_train    - Number of training set trials
%  X_test       - Test input matrix (size #cells
%                 x #time_points_per_trial * #test_set_trials)
%                 composed of horizontally concatenated population responses
%  groups_test  - Class labels of test set trials (e.g. stimulus identity);
%                 (vector of length #test_set_trials)
%  n_e_test     - Number of test set trials
% Output:
%  n_m          - Optimal number of spatiotemporal modules

n_m_max = size(X_train,1) * (size(X_train,2) / n_e_train); % Maximum number of modules

%% Find optimal number of modules
m_range = 1:n_m_max;
% Percent correct classified on the training set as a function of #modules
ctr = zeros(length(m_range),1);
% Percent correct classified on the test set as a function of #modules
cte = zeros(length(m_range),1);
for i_m = 1:length(m_range)
    n_m = m_range(i_m);
    % Decompose training set
    [Acal_train,Wst] = stnmf(X_train,n_m,n_e_train);
    % Obtain test set activation coefficients for given modules
    Acal_test = stnmf(X_test,n_m,n_e_test,Wst);
    % Process activation coefficients for classification
    predictors_train = Acal_train';
    predictors_test = Acal_test';
    [cc_train,cc_test] = ldacc(predictors_train,groups_train,predictors_test,groups_test);
    ctr(i_m) = cc_train;
    cte(i_m) = cc_test;
end

[~,i_m] = max(cte(:));
n_m = m_range(i_m);

end
