function [cur_rho, area_label] = MDSPlotInteractions(data,useProj,ndim)
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

ovr_betas = cell(N_areas, 1); %using cell arrays to help keep things straight since lots of same sizes (8 brain regions)
cur_rho = cell(N_areas, 1);
for cur_rec = 1:nrec %loop through recordings
    fprintf('\nworking on rec %d',cur_rec);
    for cur_motif = 1:nmotif %loop through motifs
        fprintf('\nworking on motif %d', cur_motif);

        %load betas
        cur_betas = loadBetas(data,cur_rec,cur_motif,ndim,area_label,useProj);

        % How many neurons per region
        N_neurons = max(cellfun(@(x) size(x, 1), cur_betas), [], 2);

        % Calculate rho's, watch out for missing data        
        for source_area = 1:N_areas,
            if isempty(cur_rho{source_area}),
                cur_rho{source_area} = NaN*ones(N_areas, N_areas, nmotif, nrec);
            end
            for pri_area = 1:N_areas,
                for sec_area = 1:N_areas,
                    % Check if either is empty, if so, then skip
                    if ~isempty(cur_betas{source_area, pri_area}) && ~isempty(cur_betas{source_area, sec_area}),
                        %temp_rho = corr(cur_betas{source_area, i}, cur_betas{source_area, j});
                        %Find the dimensions that give the maximum rhos for this pair of regions
                        %[max_rho, ~] = FindMaxRho(temp_rho, 1);

                        % Let's find the order of dimensions that gives us the highest
                        % correlation -- if dimensions are orthogonal this
                        % is equivalent to rotating
                        all_inds = perms([1:size(cur_betas{source_area, pri_area}, 2)]);
                        avg_rho = zeros(size(all_inds, 1), 1);
                        for i = 1:size(all_inds, 1),
                            temp_betas_pri = [];
                            temp_betas_sec = [];
                            for j = 1:size(cur_betas{source_area, pri_area}, 2),
                                temp_betas_pri = cat(1, temp_betas_pri, cur_betas{source_area, pri_area}(:, j));
                                temp_betas_sec = cat(1, temp_betas_sec, cur_betas{source_area, sec_area}(:, all_inds(i, j)));
                            end %dimension loop
                            
                            % Do correlation and save
                            avg_rho(i) = abs(corr(temp_betas_pri, temp_betas_sec));
                        end % index loop
                        %Take the maximum
                        [max_rho, max_rho_ind] = max(avg_rho);
                        max_rho_ind = all_inds(max_rho_ind, :);

                        % Add to our overall matrix
                        cur_rho{source_area}(pri_area, sec_area, cur_motif, cur_rec) = max_rho;
                    end
                end
            end

        end %source area loop

        % Compress all of the dimensions together
        temp_betas = cell(N_areas, 1);
        for i = 1:size(cur_betas, 1),
            for j = 1:size(cur_betas, 2),
                if (i == j), continue; end
                if size(cur_betas{i,j}, 1) == 0,
                    temp_betas{i} = cat(2, temp_betas{i}, NaN*ones(N_neurons(i), ndim));
                else
                    temp_betas{i} = cat(2, temp_betas{i}, cur_betas{i,j});
                end
            end
        end

        % Now add them to overall
        for i = 1:size(cur_betas, 1),
            if size(temp_betas{i}, 1) == 0, continue; end
            ovr_betas{i} = cat(1, ovr_betas{i}, temp_betas{i});
        end

    end %motif loop
end %rec loop

% Plot MDS network
for source_area = 1:N_areas,
    % Take the median
    temp_rho = fisherInverse(nanmedian(nanmedian(fisherZ(cur_rho{source_area}), 3), 4));
    temp_rho2 = temp_rho.^2;

    temp_rho2([1:size(temp_rho2, 1)] + size(temp_rho2, 1)*([1:size(temp_rho2, 1)]-1)) = 1;
    temp_rho2 = temp_rho2([1:(source_area-1) (source_area+1):size(temp_rho2, 1)], [1:(source_area-1) (source_area+1):size(temp_rho2, 1)]);

    PlotMDSNetwork(temp_rho2, 'CalcCorrelation', 0, 'N_areas', N_areas, 'N_dims', 1, ...
        'EdgeDataThreshold', 0.05, 'NodeEdgeColor', 'black', ...
        'Radius', 0.05, 'EdgeScaling', 10, ...
        'SourceAreaName', area_label{source_area}, 'AreaNames', area_label([1:(source_area-1) (source_area+1):end]));
end % source area loop



end %function


function [max_rho, max_rho_ind] = FindMaxRho(rho, useFisherTransform)

% Finds the maximum rhos for each pair -- does this in a holistic manner
all_inds = perms([1:size(rho, 1)]);
avg_rho = zeros(size(all_inds, 1), 1);
for i = 1:size(all_inds, 1),
    for j = 1:size(rho, 1),
        if useFisherTransform,
            if abs(rho(j, all_inds(i, j))) == 1,
                % Then we need to boost this, fisherZ doesn't give quite
                % the right answer
                avg_rho(i) = avg_rho(i) + 25;
            else
                avg_rho(i) = avg_rho(i) + fisherZ(abs(rho(j, all_inds(i, j))));
            end
        else
            avg_rho(i) = avg_rho(i) + rho(j, all_inds(i, j)).^2;
        end
    end % dimension loop
end % index loop
avg_rho = avg_rho/size(rho, 1);
if useFisherTransform,
    avg_rho = fisherInverse(avg_rho);
else
    avg_rho = sqrt(avg_rho);
end
% Find maximum rho
[max_rho, max_rho_ind] = max(avg_rho);
% Return maximum rho index
max_rho_ind = all_inds(max_rho_ind, :);

end

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




