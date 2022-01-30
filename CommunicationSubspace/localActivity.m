function [xw,yw] = localActivity(x,y,num_pcs)
%Camden - timeless
%get the weightings of each neuron that account for the most local variance
%using PCA in two populations. To just look at one population, leave y=[];
%writen this way for the communication subspace pipeline

if nargin <3; num_pcs = 1; end %number of PCs to return (def=1)

%concatentate across trials
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';

%mean substracted and pca
x = x-nanmean(x)/std(x,[],1);
xw = pca(x,'NumComponents',num_pcs);

if ~isempty(y)
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
    y = y-nanmean(y);
    yw = pca(y,'NumComponents',num_pcs);
else
    yw = [];
end

    
end %function end


