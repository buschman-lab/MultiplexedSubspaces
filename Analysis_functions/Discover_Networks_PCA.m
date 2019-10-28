function [ExpVar_all,coeff,score,mu,opts,nanpxs] = Discover_Networks_PCA(varargin)

%Set options
opts = SetAnalysisOptions();
opts = ParseOptionalInputs(opts,varargin); 

%Load the preprocessed data. Must be 1 x #rec cell arrays of data & nanpxs 
load([opts.bucket opts.base opts.data_file_name],'data','nanpxs');
fprintf('\n Calculating sPCA on block %d',opts.block);

%Get current data block
X = data{opts.block}';

%for reproducibility
rng('default'); 

%Run PCA
%Here pca is finding spatial modes, where pixels are variables and
%timepoints and observations
[coeff, score, ~, ~, ~, mu] = pca(X);

%Calculate the expvariance of each additional component
ExpVar_all = ones(1,size(score,2));
for cur_comp = 1:size(score,2)
    Xhat = repmat(mu,size(score,1),1) + score(:,1:cur_comp)*coeff(:,1:cur_comp)';    
    ExpVar_all(cur_comp) = CalculateExplainedVariance(X,(X-Xhat));
end

end %end function loop








