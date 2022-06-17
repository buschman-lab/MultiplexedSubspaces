% Does an exhaustive search of the subspace dimensions in each region that
% gives the maximum correlation across all regions. This will be the
% 'subspace network'.

%Tim Buschman - 5/11/22
% Based on code from Camden MacDowell (RunExample.m)

%% addpaths
addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));

%% load fitted models (~60 secs)
data = LoadSubspaceData('in_grouped');

% figure parameters class to make life beautiful
fp = fig_params_cortdynamics;
fprintf('Data loaded.\n')

%% Plot MDS projections of connections between nodes
useProj = 1; %def = 1; use the timepoint x trial projection along betas. 0 = just correlate betas
ndim = 5; %number of dimensions to use for exhaustive search
[rho, area_label] = MDSPlotInteractions(data,useProj,ndim);

%% Lets map the distribution of betas across different areas/subspaces
useProj = 1; %def = 1; use the timepoint x trial projection along betas. 0 = just correlate betas
ndim = 6; %number of dimensions to use for exhaustive search
[networks, area_label] = MapSubspaceNetworks(data, 'HIPP', useProj, ndim);
fprintf('\nInteractions compiled.\n')


%% Lets estimate the variance of different network structures
useProj = 1; %def = 1; use the timepoint x trial projection along betas. 0 = just correlate betas
ndim = 6; %number of dimensions to use for exhaustive search
[networks, area_label] = OptimizeNetworkVariance(data, 'HIPP', useProj, ndim);
fprintf('\nInteractions compiled.\n')



%% compile our interactions (~10 seconds) | USING A UNIFORM SAMPLE (same dimension)
useProj = 1; %def = 1; use the timepoint x trial projection along betas. 0 = just correlate betas
ndim = 4; %number of dimensions to use for exhaustive search
[rho, area_label] = CompileInteractions(data,useProj,ndim,'exhaustive');
fprintf('\nInteractions compiled.\n')

%% Now need to find the optimized networks that maximize correlation

num_areas = length(area_label);

% What is the source region?
source_area = 1;
fprintf('Creating optimized networks for subspace networks communicating with region %s (%d)\n', area_label{source_area}, source_area);

% Create an exhaustive list of subspace dimension pairings
pair_list = [];
fprintf('Creating list of networks...');
for cur_area = 1:num_areas,
    temp = [];
    for i = 0:ndim,
         temp = cat(1, temp, ones((ndim+1)^(cur_area-1), 1)*i);
    end
    pair_list = cat(2, pair_list, repmat(temp, [(ndim+1)^num_areas/length(temp) 1]));
end
fprintf('done.\n');


% What recordings and motifs do we want to look at?
selected_recordings = [];
if isempty(selected_recordings), selected_recordings = [1:size(rho,1)]; end % if no selected recordings, use them all
selected_motifs = [];
if isempty(selected_motifs), selected_motifs = [1:size(rho,2)]; end % if no selected motif, use them all
fprintf('Averaging for recordings %s and motifs %s\n', mat2str(selected_recordings), mat2str(selected_motifs));

% Loop through each one, estimate the average correlation for the network
avg_rho = NaN*ones(size(pair_list, 1), 1);
fprintf('Calculating average rhos for different networks...\n'); start_time = now;
for cur_pair = 1:size(pair_list, 1),
    % Build our list of indexes
    index_list = [];
    for cur_pair_i = 1:(num_areas-1),
        % If the pair_list includes 0 for this pair, then skip
        if pair_list(cur_pair, cur_pair_i) == 0, continue; end
        for cur_pair_j = (cur_pair_i+1):num_areas,
            % If the pair_list includes 0 for this pair, then skip
            if pair_list(cur_pair, cur_pair_j) == 0, continue; end
            % We want to grab the correlation between the specified dimensions in 
            % the pair list from region (cur_pair_i) and region (cur_pair_j)
            %
            % This should create a list with length (num_areas-1)*(num_areas-2)/2
            index_list = cat(1, index_list, ...                
                sub2ind(size(rho{1}{1}), ...
                (cur_pair_i-1)*ndim + pair_list(cur_pair, cur_pair_i), ...
                (cur_pair_j-1)*ndim + pair_list(cur_pair, cur_pair_j)));
        end 
    end

    % We're going to average over list of recordings and list of motifs
    rho_sum = 0; rho_count = 0;
    for cur_rec = 1:length(selected_recordings),
        for cur_motif = 1:length(selected_motifs),
            if isempty(index_list),
                rho_sum = 0; rho_count = 0;
            else
                %rho_sum = nansum(nansum(fisherZ(rho{selected_recordings(cur_rec), selected_motifs(cur_motif)}{source_area}(index_list))));
                rho_sum = nansum(nansum(rho{selected_recordings(cur_rec), selected_motifs(cur_motif)}{source_area}(index_list)));
                rho_count = sum(sum(~isnan(rho{selected_recordings(cur_rec), selected_motifs(cur_motif)}{source_area}(index_list))));
            end
        end % motif loop
    end %recording loop
    %avg_rho(cur_pair) = fisherInverse(rho_sum/rho_count);
    avg_rho(cur_pair) = rho_sum/rho_count;

    % Print progress
    if mod(cur_pair, floor(size(pair_list, 1)/20)) == 0,
        elapsed_time_s = (now-start_time)*24*60*60;
        fprintf('\tFinished %2.1f%% of networks...[%3.1f seconds, ~%3.1f to go]\n', cur_pair/size(pair_list, 1)*100, elapsed_time_s, elapsed_time_s/cur_pair*(size(pair_list, 1)-cur_pair));
    end
end % pair loop
fprintf('\tDone!\n');

%% Now that we've done an exhaustive search, lets iteratively find the EXCLUSIVE networks with the maximum correlation
% These DO NOT HAVE REPEATS OF DIMENSIONS
cur_avg_rho = avg_rho;
cur_pair_list = pair_list;

% Take out any NaNs
cur_pair_list = cur_pair_list(~isnan(cur_avg_rho), :);
cur_avg_rho = cur_avg_rho(~isnan(cur_avg_rho), :);

max_rho_network = []; %NaN*ones(ndim, num_areas);
max_rho_val = []; %NaN*ones(ndim, 1);
while ~isempty(cur_pair_list),
    %Find the current maximum rho
    [max_avg_rho, max_ind] = max(cur_avg_rho);

    %Save to our variables
    max_rho_val = cat(1, max_rho_val, max_avg_rho);
    max_rho_network = cat(1, max_rho_network, cur_pair_list(max_ind, :));

    max_avg_rho
    cur_pair_list(max_ind, :)

    %Now we need to remove any networks that share those dimensions
    bad_ind = zeros(length(cur_pair_list), 1);
    for cur_area = 1:num_areas,
        if cur_pair_list(max_ind, cur_area) == 0, continue; end %skip zeros
        bad_ind = or(bad_ind, cur_pair_list(:, cur_area) == cur_pair_list(max_ind, cur_area));
    end
    cur_avg_rho = cur_avg_rho(~bad_ind);
    cur_pair_list = cur_pair_list(~bad_ind, :);
end



