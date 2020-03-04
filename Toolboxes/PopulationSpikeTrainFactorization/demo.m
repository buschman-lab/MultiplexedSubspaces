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

% This script demonstrates how to apply the Population Spike Train
% Factorization Toolbox. It loads an example dataset, converts its spike
% trains spike count matrices and applies spatiotemporal NMF and
% space-by-time NMF to evaluate spatiotemporal and space-by-time decoding
% performances.


% Clear all variables
clear all;
% Load example dataset
disp('Loading data...');
load('example_data.mat');
% Bin size in milliseconds
bin_size = 10;
% Length of a trial in milliseconds
len_trial = 300;

%% Transform spike trains to spike count matrices
disp('Transforming spike trains to count matrices...');
counts = cell(n_stimuli,1);
for i = 1:n_stimuli
    counts{i} = count_spikes(Spikes{i},bin_size,len_trial);
end

clear i;

%% Randomly separate data into training set and test set
p = randperm(n_trials); % Random permutation of trials
% Training indices
ind_train = p(1:ceil(n_trials/2));
% Test indices
ind_test = p((ceil(n_trials/2)+1):end);

% Total number of training samples
n_e_train = n_stimuli*length(ind_train);
% Total number of test samples
n_e_test = n_stimuli*length(ind_test);
% Number of bins per trial
n_bins = ceil(len_trial/bin_size);

clear p;

%% Build overall data matrices
disp('Concatenating data...');
% Training set
X_train = zeros(n_neurons,n_e_train*n_bins);
offset = 0;
for i = 1:n_stimuli
    for j = 1:length(ind_train)
        X_train(:,offset+(1:n_bins)) = counts{i}{ind_train(j)};
        offset = offset + n_bins;
    end
end
% Training class labels
groups_train = ceil((1:n_e_train)' / length(ind_train));

% Test set
X_test = zeros(n_neurons,n_e_test*n_bins);
offset = 0;
for i = 1:n_stimuli
    for j = 1:length(ind_test)
        X_test(:,offset+(1:n_bins)) = counts{i}{ind_test(j)};
        offset = offset + n_bins;
    end
end
% Test class labels
groups_test = ceil((1:n_e_test)' / length(ind_train));

clear offset i j;

%% Apply spatiotemporal NMF to factorize data
disp(['Optimizing number of spatiotemporal modules.'...
    ' This can take a while...']);
% Find optimal number of spatiotemporal modules
n_m = select_n_m(X_train,groups_train,n_e_train,X_test,groups_test,...
    n_e_test);
% Obtain spatiotemporal modules from training set
[Acal_train,Wst] = stnmf(X_train,n_m,n_e_train);
disp(['Found ' int2str(n_m) ' spatiotemporal modules.']);
% Obtain activation coefficients from test set for given modules
Acal_test = stnmf(X_test,n_m,n_e_test,Wst);
% Process activation coefficients for classification
predictors_train = Acal_train';
predictors_test = Acal_test';
% Get classification performance on training and test sets
[cc_st_train,cc_st_test] = ldacc(predictors_train,groups_train,...
    predictors_test,groups_test);

clear n_m Acal_train Acal_test predictors_train predictors_test;

%% Apply space-by-time NMF to factorize data
disp(['Optimizing numbers of temporal and spatial modules.'...
    ' This can take a while...']);
% Find optimal numbers of temporal and spatial modules
[n_tm,n_sm] = select_n_tm_n_sm(X_train,groups_train,n_e_train,X_test,...
    groups_test,n_e_test);
% Obtain temporal and spatial modules from training set
[Acal_train,Wi,Wb] = sbtnmf(X_train,n_tm,n_sm,n_e_train);
disp(['Found ' int2str(n_tm) ' temporal and ' int2str(n_tm)...
    ' spatial modules.']);
% Obtain activation coefficients from test set for given modules
Acal_test = sbtnmf(X_test,n_tm,n_sm,n_e_test,Wi,Wb);
% Process activation coefficients for classification
predictors_train = zeros(n_e_train,n_tm*n_sm);
for i = 1:n_e_train
    predictors_train(i,:) = reshape(Acal_train(:,:,i),1,n_tm*n_sm);
end
predictors_test = zeros(n_e_test,n_tm*n_sm);
for i = 1:n_e_test
    predictors_test(i,:) = reshape(Acal_test(:,:,i),1,n_tm*n_sm);
end
% Get classification performance on training and test sets
[cc_sbt_train,cc_sbt_test] = ldacc(predictors_train,groups_train,...
    predictors_test,groups_test);

clear n_sm n_tm Acal_train Acal_test predictors_train predictors_test i;

%% Show results
disp([]);
disp(['Spatiotemporal classification performance on training set: '...
    num2str(cc_st_train) '%']);
disp(['Spatiotemporal classification performance on test set: '...
    num2str(cc_st_test) '%']);
disp(['Space-by-time classification performance on training set: '...
    num2str(cc_sbt_train) '%']);
disp(['Space-by-time classification performance on test set: '...
    num2str(cc_sbt_test) '%']);

