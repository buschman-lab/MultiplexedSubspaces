function [networks, area_label] = OptimizeNetworkVariance(data, source_area_str, useProj, ndim)
%Camden MacDowell - timeless
%NOTE; beta direction is random, so this returns absolute value of
%correlations

if nargin <2; useProj = 1; end %def = 1; use the timepoint/trial projection along each beta (slower). 0 = just correlate betas
if nargin <3; ndim = 4; end %number of dimensions to use

nrec = numel(data); %num  recordings
nmotif = size(data{1},2); %num motifs;

%full list of areas (recs 1,2,5,6 have all areas; 3 is missing RSP and BFD. 4 if missing RSP)
area_label = data{1}(1).area_label;
N_areas = length(area_label);
cur_source_area = find(strcmpi(area_label, source_area_str), 1, 'first'); %index of current source area
if isempty(cur_source_area), error('Bad source area passed.'); end

networks = cell(nrec,nmotif); %using cell arrays to help keep things straight since lots of same sizes (8 brain regions)
for cur_rec = 1:nrec %loop through recordings
    fprintf('working on rec %d\n',cur_rec);
    for cur_motif = 1:nmotif %loop through motifs
        fprintf('\tworking on motif %d\n',cur_motif);
        %load betas
        betas = loadBetas(data,cur_rec,cur_motif,ndim,area_label,useProj);
        fprintf('\t\tBetas loaded.\n');

        % Grab the betas from the source area alone
        b = betas(cur_source_area, :);
        
        % What dimensions should we use?
        dim_list = ones(N_areas, 1);
        dim_list(cur_source_area) = 0;
        cur_betas = [];
        for dim_ind = 1:length(dim_list),
            if dim_list(dim_ind),
                cur_betas = cat(2, cur_betas, b{dim_ind});
            end
        end

        % Now do SVD to get all of the dimensions together
        [U, S, V] = svd((1/(size(cur_betas, 1) - 1))*(cur_betas'*cur_betas));

        % Plot the distribution of PEV of each dimension
        figure; plot(cumsum(diag(S).^2./(sum(diag(S).^2)))); hold on;
        bar(diag(S).^2./(sum(diag(S).^2)));

        figure;
        for i = 1:(N_areas-1),
            plot([1:ndim] + (i-1)*(ndim+1), abs(U((i-1)*ndim + [1:ndim], 1:4))); 
            hold on;
        end

    end %motif loop
end %rec loop

end %function

function PlotPairwiseLineFits(b),

dim_list = [NaN 1 1 1 1 1 1 1];
cur_betas = [];
for dim_ind = 1:length(dim_list),
    if ~isnan(dim_list(dim_ind)),
        cur_betas = cat(2, cur_betas, b{dim_ind}(:, dim_list(dim_ind)));
    end
end

% Subtract the mean from everything to get correlations/fits to behave
% nicely
beta_offset = nanmean(cur_betas, 1);
cur_betas = cur_betas - repmat(beta_offset, [size(cur_betas, 1) 1]);

figure;
for i = 1:size(cur_betas, 2),
    for j = 1:size(cur_betas, 2),
        subplot(size(cur_betas, 2), size(cur_betas, 2), (i-1)*size(cur_betas, 2) + j);
        plot(cur_betas(:, i), cur_betas(:, j), '.');

        % Calculate the correlation
        cur_rho = corr(cur_betas(:, i), cur_betas(:, j), 'type', 'pearson');

        % Let's fint a line of maximal explained variance -- test code for
        % ND in future
        cur_beta_mat = cur_betas(:, [i j]);
        [U, S, V] = svd((1/size(cur_beta_mat, 1))*(cur_beta_mat'*cur_beta_mat));
        % The first dimension should be our line of best fit
        beta_line_fit = U(:, 1);
        %Plot it
        v = axis;
        hold all; plot([v(1) v(2)], beta_line_fit(2)/beta_line_fit(1)*[v(1) v(2)], 'r-');
        %And the other dimensions should be the residuals from that fit
        cur_resid = U(:, 2:end)'*cur_beta_mat';
        SST = sum(sum(cur_beta_mat.^2));
        SSE = sum(sum(cur_resid.^2));
        est_r2 = 1 - SSE/SST;

        % Use PCA as a gut check
        [coeff,score,latent,tsquared,explained,mu] = pca(cur_beta_mat);
        % Take the first PC and use it as our line of best fit
        pc_line_fit = coeff(:, 1);
        hold all; plot([v(1) v(2)], pc_line_fit(2)/pc_line_fit(1)*[v(1) v(2)], 'g-');
        %And the other dimensions should be the residuals from that fit
        cur_resid = coeff(:, 2:end)'*cur_beta_mat';
        SST = sum(sum(cur_beta_mat.^2));
        SSE = sum(sum(cur_resid.^2));
        est_pca_r2 = 1 - SSE/SST;


        % Add title
        title(gca, sprintf('rho = %3.2f, r^2 = %3.2f (%3.2f)', cur_rho, est_r2, est_pca_r2));
    end
end
end %Pairwise line fits

%% Function to load betas

function [betas,targidx] = loadBetas(data,cur_rec,cur_motif,ndim,area_label,useProj)
%retruns betas which is an 8x8 cell array. Each cell cointains matrix
%of betas weights of a region (nuerons x dimensions)
% useProj, then returns (timepoints x trials) x dimensions matrix of
% activtiy in a region projected along those beta weights
%rows of beta are predictor regions. columns are targets
%i.e., row 1 of beta column 2 is (neurons x dimensions) betas for Hipp predicting area MOs
%Gutchecks:
%betas should have nothing on the diagnol
%betas for rec three rows and columns 4/6 should be empty and rec 4 row/column 3

b = data{cur_rec}(cur_motif).rrr_B; %each cell is beta weights of n-1 regions predicting each target region
cur_areas = data{cur_rec}(cur_motif).area_label; %the list of target regions

%just keep desired dimensions
b = cellfun(@(x) x(:,1:ndim),b,'UniformOutput',0);

%reorganize betas so each cell is the betas for a given area when predicting all other areas
g = data{cur_rec}(cur_motif).grouping; %each cell is region of beta weights of n-1 regions predicting each target region
betas = cell(numel(area_label),numel(area_label));
targidx = ismember(area_label,cur_areas); %get true index for targets (out of 8 areas, to match across recs)
for cur_a = 1:numel(cur_areas)
    idx = strcmp(area_label,cur_areas(cur_a)); %get true index (out of 8 areas, to match across recs)
    if ~isempty(idx)
        betas(idx,targidx) = arrayfun(@(n) b{n}(g{n}==cur_a,:),1:numel(g),'UniformOutput',0);
    end
end

%if projecting
if useProj == 1
    y = data{cur_rec}(cur_motif).area_val;
    y = cellfun(@(x) normalizeToBaseline(x,[1:2],'mean'), y,'UniformOutput',0);%normalize to baseline
    y = cellfun(@(x) x(:,3:end,:), y,'UniformOutput',0);%use non-baseline
    y = cellfun(@(x) x-nanmean(x,3), y,'UniformOutput',0);%subtract the psth
    y = cellfun(@(x) reshape(x,[size(x,1),size(x,2)*size(x,3)])', y,'UniformOutput',0);%reshape into trial/timepoint x neuron

    %reshape to match predictors in beta
    yy = cell(8,1);
    yy(targidx)=y;
    yy = repmat(yy,1,8);
    %temporily replace emptys with nan for multiplication
    badidx = cellfun(@isempty,betas);
    yy(badidx) = {NaN};
    betas(badidx) = {NaN};
    betas = cellfun(@(x,y) x*y, yy,betas,'UniformOutput',0);
    betas(badidx) = {[]};

end %proj if

end









