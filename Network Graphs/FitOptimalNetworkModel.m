function [proj_data, proj_resid, proj_vect, ovr_SSE, ovr_SST, ovr_model_comp_stat, dims] = FitOptimalNetworkModel(cur_data, model_comp_stat, disp)
% Find the variables that maximizes the explained variance across all of
% our data points. Everything will be assumed to have a core linear dimension (1D)
% because that is what we are fitting from the ridge regression.

if nargin<2, model_comp_stat = 'BIC'; end
if nargin<3, disp = 1; end

% What is the tolerance in the best model that we are okay with
model_comp_stat_stop_tol = 100*eps; %a bit above machine precision for now
max_iters = 10000; %maximum number of iterations

% Display figure
if disp, figure; end

%Check to see if we have any NaNs
if any(isnan(cur_data(:))), error('Data must not contain NaNs. Don''t know what to do with those...'); end

% Zero-mean our data
cur_data = cur_data - repmat(mean(cur_data, 1), [size(cur_data, 1) 1]);

N_dim = size(cur_data, 2);
ovr_N = numel(cur_data);

% We are going to start with every column of our matrix being independent
dims = mat2cell([1:N_dim], 1, ones(N_dim, 1));

% Initialize with our null model
proj_data = cur_data;
cur_model_comp_stat = Inf;
iter_count = 0;
prev_model_comp_stat = Inf;
ovr_model_comp_stat = [];
ovr_SSE = [];
ovr_SST = [];

% Do this repeatedly, as long as it improves our model
while (iter_count <= max_iters) && (((prev_model_comp_stat - cur_model_comp_stat) >= model_comp_stat_stop_tol) || (iter_count == 0)),
    % Update our previous model comparison stat
    prev_model_comp_stat = cur_model_comp_stat;

    % Combine dims of most correlated dimensions
    dims = CombineDims(proj_data, dims);

    % Do it all again!
    [proj_data, proj_resid, proj_vect, SSE, SST] = ProjectData(cur_data, dims);
    % Get the overall number of parameters
    k = 0;
    for i = 1:length(proj_vect),
        k = k + numel(proj_vect{i}) - 1;
    end
    ovr_SSE = cat(1, ovr_SSE, SSE);
    ovr_SST = cat(1, ovr_SST, SST);

    % Estimate model comparison statistics
    cur_model_comp_stat = ModelCompStat(N_dim, k, SST, SSE, model_comp_stat);
    iter_count = iter_count + 1;
    ovr_model_comp_stat = cat(1, ovr_model_comp_stat, [iter_count cur_model_comp_stat]);

    if disp,
        subplot(1,2,1); cla;
        disp_mat = zeros(length(dims), size(cur_data, 2));
        for i = 1:length(dims),
            disp_mat(i, dims{i}) = 1;
        end
        imagesc(disp_mat);
        title(sprintf('Dimension Map', model_comp_stat));

        subplot(1,2,2); cla;
        plot(ovr_model_comp_stat(2:end, 1), ovr_model_comp_stat(2:end, 2), 'k-');
        xlabel('# of Iterations'); ylabel(model_comp_stat);
        title(sprintf('%s Over Iterations', model_comp_stat));
        drawnow;
    end
end % loop

end 

%% Function to combine two most correlated dimensions
function dims = CombineDims(proj_data, dims),
    % Find our most correlated dimensions
    rho = corr(proj_data);
    rho = triu(rho.^2, 1); 
    [~, max_ind] = max(rho(:));
    [dim_x, dim_y] = ind2sub(size(rho), max_ind);

    % Combine these and get rid of the old dimension
    dims{dim_x} = cat(2, dims{dim_x}, dims{dim_y});
    dims{dim_y} = [];
    dims = dims(~cellfun(@isempty, dims));
end

%% Function to do projections
function [proj_data, proj_resid, proj_vect, SSE, SST] = ProjectData(cur_data, dims),
    % For each index in the dims variable, we are going to take the line of
    % maximal variance

    %Initialize output variables
    proj_data = NaN*ones(size(cur_data, 1), length(dims));
    proj_resid = cell(length(dims), 1);
    proj_vect = cell(length(dims), 1);
    SSE = 0; SST = 0;

    % Loop through each dimension, fitting and adding to model goodness of fit
    for cur_dim = 1:length(dims),
        % Take just the data from this dimension
        dim_data = cur_data(:, dims{cur_dim});

        % Fit a line that maximizes variance explained
        dim_covar_mat = (1./(size(dim_data, 1) - 1))*(dim_data'*dim_data);
        [U, ~, ~] = svd(dim_covar_mat);
        proj_vect{cur_dim} = U(:, 1);

        % Project onto the fit vector
        proj_data(:, cur_dim) = dim_data*U(:, 1);

        % Get the residuals
        proj_resid{cur_dim} = fit_residuals(dim_data, U(:, 1));

        % Estimate errors
        if size(dim_data, 2) == 1,
            SSE = SSE + sum(sum((dim_data - repmat(mean(dim_data, 1), [size(dim_data, 1) 1])).^2, 1), 2);
        else
            SSE = SSE + sum(sum((proj_resid{cur_dim} - repmat(mean(proj_resid{cur_dim}, 1), [size(proj_resid{cur_dim}, 1) 1])).^2, 1), 2);
        end
        SST = SST + sum(sum((dim_data - repmat(mean(dim_data, 1), [size(dim_data, 1) 1])).^2, 1), 2);

    end %dim loop

end % ProjectData function


%% Function to remove the contribution of the current vector
function cur_data = fit_residuals(cur_data, z),
    znorm=z/norm(z);
    cur_data = cur_data - (cur_data*znorm)*znorm';

%      z=cur_data*z;
%      znorm=z/norm(z);
%            
%      % perform deflation
%      cur_data= cur_data - znorm*(znorm'*cur_data);
end %deflate_data function

%% Function to estimate model comparison statistics
function start_model_comp_stat = ModelCompStat(n, k, SSE, SST, model_comp_stat),
    % Pass:
    %   N = total number of data points
    %   k = number of parameters
    %   SSE = sum squared error of fit
    %   SST = sum squared total error
    %   model_comp_state = what model comparison statistic to use

    switch model_comp_stat
        case 'MSE'
            start_model_comp_stat = -SSE/n;
        case 'r2'
            start_model_comp_stat = 1 - SSE./SST;
        case 'AIC'
            start_model_comp_stat = n*log(SSE/n) + 2*k;
        case 'AICc'
            start_model_comp_stat = n*log(SSE/n) + 2*k + (2*k^2 + 2*k)/(n-k-1);
        case 'BIC'
            start_model_comp_stat = n*log(SSE/n) + k*log(n);
        case 'adjR2'
            start_model_comp_stat = -(1 - ((n-1)/(n-k))*(1 - (1 - SSE/SST))); %CHECK THIS!! -- this should actually be 1-adjR2 to match other statistics in that minimizing should be good
        case 'Amemiya'
            start_model_comp_stat = ((n+k)/(n-k))*(SSE/n);
        otherwise
            error('Unknown model comparison statistic passed.');
    end

end 

