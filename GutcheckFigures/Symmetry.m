function Symmetry(stack,opts)

if nargin <3
    % load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\ProcessedData_Hemo\432-10-18-2019_1dff_combined.mat','stack','opts');    
end

% just take a chunk: 
dff = stack(:,:,1:opts.fps*5*60);

%spatially bin to 68x68
dff = SpatialBin(dff,2,[],1);

% filter
dff_filt = filterstack(dff, 13, [0.1 4], 'lowpass', 1, 0);

%split into right and left hemispheres and flip to match
rh_dff = dff_filt(:,2:34,:);
lh_dff = fliplr(dff_filt(:,36:end,:));

%double check to make sure we okay
% imagesc(nanmean(rh_dff,3))
% imagesc(nanmean(lh_dff,3))
% imshowpair(mean(rh_dff,3),mean(lh_dff,3))

[X,Y,Z] = size(rh_dff);
rh_dff = reshape(rh_dff,X*Y,Z);
lh_dff = reshape(lh_dff,X*Y,Z);
bad_pxl = nanvar(lh_dff,[],2) <= eps | nanvar(rh_dff,[],2)<=eps;
rh_dff(bad_pxl,:) = [];
lh_dff(bad_pxl,:) = [];

N = size(rh_dff,1);

rho = NaN(N,1);
d = NaN(N,1);

for i = 1:N
    rho(i) = corr(rh_dff(i,:)',lh_dff(i,:)');
    d(i) = sqrt(nanmean((rh_dff(i,:)-lh_dff(i,:)).^2));
end


rho_recon = NaN(1,numel(bad_pxl));
rho_recon(bad_pxl==0) = rho;
rho_recon = reshape(rho_recon,[X,Y,1]);
figure; hold on; 
imagesc(rho_recon,[-1 1]); axis equal; colormap magma
set(gca,'ydir','reverse','box','off','xlim',[1,Y],'ylim',[1,X],'xtick',[],'ytick',[])
colorbar

d_recon = NaN(1,numel(bad_pxl));
d_recon(bad_pxl==0) = d;    
d_recon = reshape(d_recon,[X,Y,1]);
figure; hold on; 
imagesc(d_recon); axis equal
set(gca,'ydir','reverse','box','off','xlim',[1,Y],'ylim',[1,X],'xtick',[],'ytick',[])
colorbar

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




