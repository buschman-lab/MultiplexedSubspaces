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

function [n_tm,n_sm] = select_n_tm_n_sm(X_train,groups_train,n_e_train,X_test,groups_test,n_e_test)
% Selects the optimal numbers of temporal and spatial modules.
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
%  n_tm         - Optimal number of temporal modules
%  n_sm         - Optimal number of spatial modules

n_tm_max = size(X_train,2) / n_e_train; % Maximum number of temporal modules
n_sm_max = size(X_train,1); % Maximum number of spatial modules

%% Find optimal number of modules
sm_range = 1:n_sm_max;
tm_range = 1:n_tm_max;
% Percent correct classified on the training set as a function of #modules
ctr = zeros(length(sm_range),length(tm_range));
% Percent correct classified on the test set as a function of #modules
cte = zeros(length(sm_range),length(tm_range));
for i_sm = 1:length(sm_range)
    n_sm = sm_range(i_sm); % Number of spatial modules
    for i_tm = 1:length(tm_range)
        n_tm = tm_range(i_tm); % Number of temporal modules

        % Decompose training set
        [Acal_train,Wi,Wb] = sbtnmf(X_train,n_tm,n_sm,n_e_train);
        % Obtain test set activation coefficients for given modules
        Acal_test = sbtnmf(X_test,n_tm,n_sm,n_e_test,Wi,Wb);
        % Process activation coefficients for classification
        predictors_train = zeros(n_e_train,n_tm*n_sm);
        for i_s = 1:n_e_train
            predictors_train(i_s,:) = reshape(Acal_train(:,:,i_s),1,n_tm*n_sm);
        end
        predictors_test = zeros(n_e_test,n_tm*n_sm);
        for i_s = 1:n_e_test
            predictors_test(i_s,:) = reshape(Acal_test(:,:,i_s),1,n_tm*n_sm);
        end
        [cc_train,cc_test] = ldacc(predictors_train,groups_train,predictors_test,groups_test);
        ctr(i_sm,i_tm) = cc_train;
        cte(i_sm,i_tm) = cc_test;
    end
end

[~,i] = max(cte(:));
[i_sm,i_tm] = ind2sub(size(cte),i);
n_tm = tm_range(i_tm);
n_sm = sm_range(i_sm);

end

