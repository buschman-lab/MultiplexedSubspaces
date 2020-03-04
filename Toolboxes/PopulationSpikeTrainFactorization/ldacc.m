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

function [cc,cc_test] = ldacc(predictors,groups,predictors_test,groups_test)
% Classifies samples by means of linear discriminant analysis and calculates
% percent correct classified.
% Input arguments:
%  predictors      - Predictors of each trial (e.g. population responses)
%                    (size #trials x #predictors_per_trial)
%  groups          - Class labels of each trial (e.g. stimulus identity);
%                    (vector of length #trials)
% Optional arguments:
%  predictors_test - Predictors of test set trials (e.g. population responses)
%                    (size #test_set_trials x #predictors_per_trial)
%  groups_test     - Class labels of test set trials (e.g. stimulus identity);
%                    (vector of length #test_set_trials)
% Output:
%  cc              - Percentage correct classified for predictors and groups
%  cc_test         - Percentage correct classified for the test set using
%                    the LDA training from the predictors and groups; empty
%                    if predictors_test and groups_test are not given

%% Apply linear discriminant analysis
% Exclude zero variance predictors
nz = var(predictors) > eps;
if (~any(nz))
    error('ldacc: No informative predictor found');
end

predictors = predictors(:,nz);
cls = ClassificationDiscriminant.fit(predictors,groups,'discrimType','pseudoLinear');
result_training = predict(cls,predictors);
cc = sum(result_training == groups) ./ length(groups) * 100;

%% Eventually evaluate test set performance
if (nargin > 2)
    % Apply linear discriminant analysis
    predictors_test = predictors_test(:,nz);
    result_test = predict(cls,predictors_test);
    cc_test = sum(result_test == groups_test) ./ length(groups_test) * 100;
else
    cc_test = [];
end

end

