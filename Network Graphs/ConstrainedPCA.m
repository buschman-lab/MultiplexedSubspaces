function [cur_vect, exp_var] = ConstrainedPCA(cur_data, constraints, lambda, norm_type, disp)
% Find a dimension of maximal explained variance with constraints on
% weights

if nargin<3, lambda = 10; end
if nargin<4, norm_type = 'L1'; end
if nargin<5, disp = 1; end

% Display figure
if disp, figure; end

%Constraints -- should be a cell array with each cell containing a list of
%dimensions. Only one of the dimensions within this list will be used at
%any point in time (others will be zero).
%constraints = {[1:6], [7:12], [13:18], [19:24], [25:30], [31:36], [37:32]};

% Zero-mean our data
cur_data = cur_data - repmat(nanmean(cur_data, 1), [size(cur_data, 1) 1]);

% Start with the first PC as our initial guess
covar_mat = (1./(size(cur_data, 2) - 1)) * (cur_data'*cur_data);
[U, S, V] = svd(covar_mat);
cur_vect = U(:, 1);
if disp, plot(cur_vect, 'r-'); hold all; end

if strcmpi(norm_type, 'GreedySVD'),
    % Sequentially drop the least important component from the data and
    % re-fit with the remaining data
    
    %Start with everybody being in our SVD
    good_ind = ones(size(cur_data, 2), 1);
    good_constraints = cell(size(constraints));
    for i = 1:length(constraints),
        good_constraints{i} = ones(length(constraints{i}), 1);
    end

    stop = 0;
    while ~stop,
        stop = 1;
        
        %Loop through constraint groups, dropping minimum each time
        for i = 1:length(constraints),
            % Find the minimum contributor
            

            %If we have more than one factor left for this constraint group, then continue
            if sum(good_constraints{i}) > 1,
                stop = 0;
            end
            
        end %constraint group loop

        

    end

else
    % Use a least-squares fit to estimate with some costs

    %Zero out all of the non-maximum within our constraints
    for i = 1:length(constraints),
        [~, max_ind] = max(abs(cur_vect(constraints{i})));
        cur_vect(constraints{i}(constraints{i} ~= constraints{i}(max_ind))) = 0;
    end
    if disp, plot(cur_vect, 'b-'); hold all; end

    %Create function to estimate SSE
    SSE = @(x) (sum(sum(deflate_data(cur_data, x).^2)));

    % Create function for general norm
    general_norm = @(x) constraint_norm(x, {[1:length(x)]}, lambda, norm_type);

    %Create optimization function -- we want to MINIMIZE this
    optim_func = @(x) (SSE(x) + constraint_norm(x, constraints, lambda, norm_type) + general_norm(x));
    %optim_func = @(x) SSE(x);

    % Optimize the function
    optim_opts = optimset('Display','final', 'MaxFunEvals', 10^10, 'MaxIter', 10^10, 'TolX', 10^-6, 'TolFun', 10^-6);
    %[cur_vect, fval, exitflag, output] = fminsearch(optim_func, cur_vect, optim_opts);
    cur_vect = lsqnonlin(optim_func, cur_vect, -ones(size(cur_vect)), ones(size(cur_vect)), optim_opts);
    cur_vect = cur_vect./norm(cur_vect, 2);
    if disp, plot(cur_vect, 'g-'); hold all; end

end

% Calculate the explained variance and return
SST = sum(sum(cur_data.^2));
SSE = SSE(cur_vect);
exp_var = 1 - SSE/SST;

end 

%% Function to remove the contribution of the current vector
function cur_data = deflate_data(cur_data, z)
%z=cur_data*cur_vect;
znorm=z/norm(z);

cur_data = cur_data - (cur_data*znorm)*znorm';
end %deflate_data function

%% Function to estimate constraint norm
function cost = constraint_norm(cur_vect, constraints, lambda, norm_type),

cost = 0;
switch norm_type
    case 'L0'
        for i = 1:length(constraints),
            cost = cost + lambda*sum(abs(cur_vect(constraints{i})) > 0);
        end
    case 'L1'
        for i = 1:length(constraints),
            cost = cost + lambda*sum(abs(cur_vect(constraints{i})));
        end
    case 'L2'
        for i = 1:length(constraints),
            cost = cost + lambda*sum(cur_vect.^2);
        end
    case 'NonMaxL1'
         for i = 1:length(constraints),
            temp_cost = abs(cur_vect(constraints{i}));
            cost = cost + lambda*(sum(temp_cost) - max(temp_cost));
         end
     case 'NonMaxL2'
         for i = 1:length(constraints),
            temp_cost = cur_vect(constraints{i}).^2;
            cost = cost + lambda*(sum(temp_cost) - max(temp_cost));
        end
    otherwise
        error('Unknown cost function.')
end %switch
end % constraint_norm function