function IndividualMotifsPerCluster(data_dir,basis_name,save_dir,convolv_flag,W)
%camden macdowell - timeless

% if nargin <1; data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution'; end
% if nargin <2; basis_name = 'ALL'; end
% if nargin <3 save_dir = [data_dir filesep 'GutCheckData']; end
if ~exist(save_dir)
    mkdir(save_dir)
end

%load the cluster
load([data_dir filesep basis_name filesep 'basismotifs.mat'],'W_basis','noise_clusters','core_comm_idx');

if isempty(W)
    %load all the individual motifs
    file_list = GrabFiles('\w*chunk\w*',0,{data_dir});
    gp = general_params_vpa;
    fp = fig_params;
    %load data
    temp = cellfun(@(x) load(x,'W','nanpxs'),file_list,'UniformOutput',0);
    W = cellfun(@(x) x.W, temp,'UniformOutput',0);
    nanpxs = cellfun(@(x) x.nanpxs, temp,'UniformOutput',0);

    %Recondition W
    for i = 1:numel(W)    
       temp = zeros(gp.pixel_dim(1)*gp.pixel_dim(1),size(W{i},2),size(W{i},3));
       idx = ones(size(temp,1),1);
       idx(nanpxs{i})=0; %NOTE, NAN Will only use the same across all pixels, zeros will smooth essentailly for clustering 
       temp(idx==1,:,:) = W{i};
       W{i} = temp;
    end

    W = cat(2,W{:});
    nanpxs = find(nanvar(reshape(W,[size(W,1),size(W,2)*size(W,3)]),[],2)<=eps);
    W(nanpxs,:,:) = [];

    if iscell(W) %combine mutliple CNMF Fits stored in a cell array
       W = cat(2,W{:});
    end

    %Remove empty motifs 
    W = RemoveEmptyMotifs(W);
    W = MaskTensor(W,nanpxs,[gp.originaldimensions(1)*gp.originaldimensions(2),size(W,2),size(W,3)]); 
end

%create example gifs
[P,K,L] = size(W_basis);
for i = 1:K      
   if convolv_flag
       mu = 3;
       sd = 0.5;x = linspace(1,5,30);
       y = 1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2));
       temp = tensor_convolve(W_basis(:,i,:),cat(2,y,zeros(1,15)));
   else
       temp = squeeze(W_basis(:,i,:));
   end
   
   temp = reshape(temp,[68 68 size(temp,2)]);
   if ismember(i,noise_clusters)
      MotifToGif(temp,[save_dir filesep sprintf('BM_noise%d.gif',i)],'limit',100); 
   else
      MotifToGif(temp,[save_dir filesep sprintf('BM%d.gif',i)],'limit',100); 
   end
   
%    Get up to 20 random motifs that are part of this cluster
   temp_idx = core_comm_idx{i};
   if numel(temp_idx)>5
       temp_idx = temp_idx(randperm(numel(temp_idx),5));
   end
   indi_motif = W(:,temp_idx,:);
   for j = 1:size(indi_motif,2)
       temp = squeeze(indi_motif(:,j,:));
       temp = reshape(temp,[68 68 size(temp,2)]);
       if ismember(i,noise_clusters)
           MotifToGif(temp,[save_dir filesep sprintf('BM_noise%d_Example%d.gif',i,j)],'limit',100); 
       else
           MotifToGif(temp,[save_dir filesep sprintf('BM%d_Example%d.gif',i,j)],'limit',100); 
       end
   end   
end %basis loop
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

