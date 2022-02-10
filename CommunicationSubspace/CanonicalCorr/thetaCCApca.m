function thetaCCApca(x,x_coef,num_cvs, num_pcs)
%Camden - timeless
%x is cell array of the neuron x sig CVs from CCA
%x_coef is cel array the neuron x PCs from PCA

if nargin <3; num_cvs = 1; end
if nargin <4; num_pcs = 1; end

%future camden note; it would be weird to use multiple PCS.... since PCS
%are orthogonal.... need to think more on this. 
theta = NaN(1,numel(x)); 
for i = 1:numel(x)
    if ~isempty(x{i}) %if a sig CV was found
       theta(i) = AngleBetweenWeights(x{i}(:,1:num_cvs),x_coef{i}(:,1:num_pcs),'none');
    end
end


end %function end




