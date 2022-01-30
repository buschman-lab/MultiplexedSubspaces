function [stat,stat_mat,g] = CompareSubspacesWithinAPopulation(x,grp,type)
%Camden MacDowell - timeless
%Get the angles between unique subspaces within a brain area. grp is the
%brain area label
%type is a string specifying the comparison to make

g = unique(grp);

stat = cell(1,numel(g));
stat_mat = cell(1,numel(g));  %also return as a matrix for each area
for i = 1:numel(g)
    pair_idx = nchoosek(find(grp==g(i)),2);
    %loop through each pair
    for j = 1:size(pair_idx,1)
        switch type
            case 'full_angle' %uses all CVs
                stat{i}(j) = AngleBetweenWeights(x{pair_idx(j,1)},x{pair_idx(j,2)},'none'); 
                stat_mat{i}(pair_idx(j,1),pair_idx(j,2)) = stat{i}(j);
            case 'angle' %uses just the first CV
                stat{i}(j) = AngleBetweenWeights(x{pair_idx(j,1)}(:,1),x{pair_idx(j,2)}(:,1),'none'); 
                stat_mat{i}(pair_idx(j,1),pair_idx(j,2)) = stat{i}(j);
            case 'rho' %uses all CVs
                stat{i}(j) = fisherZ(corr(x{pair_idx(j,1)}(:,1),x{pair_idx(j,2)}(:,1)));
                stat_mat{i}(pair_idx(j,1),pair_idx(j,2)) = stat{i}(j);
            otherwise 
                error('unknown type');
        end %switch
    end %pair
end %g


end %function