function [rho_in, rho_out] = CraniotomyContraCorrelation(dff,probe_coords)
%Camden MacDowell - timeless
%returns the correlation between pixels within craniotomy to contralateral
%hemisphere vs other ipsilateral pixels vs contra hemi. 
%probe_coords is cell array of electrode insertion points
%dff can be any 68 x 68 x time pixel trace

%mask is not perfectly centered (which is deliberate for earlier preprocessing legacy), but it means we have to shift it on this end.
R_hemi = NaN(size(dff,1),size(dff,2),size(dff,3)); %convert back to full size
R_hemi(:,37:end,:) = fliplr(dff(:,3:34,:));
R_hemi = reshape(R_hemi,[size(R_hemi,1)^2,size(R_hemi,3)]);
L_hemi = NaN(size(dff,1),size(dff,2),size(dff,3)); %convert back to full size
L_hemi(:,37:end,:) = dff(:,37:end,:);
L_hemi = reshape(L_hemi,[size(L_hemi,1)^2,size(L_hemi,3)]);

mask = zeros(size(dff,1),size(dff,2));

within_mask = cell(1,4);
for i = 1:numel(probe_coords)
    temp = probe_coords{i};
    within_mask{i} = mask;
    within_mask{i}(temp(1,2)+[-2:2],temp(1,1)+[-2:2])=1; %remember, x and y need to be flipped from what you would expect
    within_mask{i} = reshape(within_mask{i},[size(mask,1)^2,size(mask,3)]); %flatten
end
y = size(L_hemi,1);
%contralateral correlation per pixel
rho = NaN(y,1);
for i = 1:y
    rho(i) = corr((R_hemi(i,:)-nanmean(R_hemi(i,:)))',(L_hemi(i,:)-nanmean(L_hemi(i,:)))');
end
rho = fisherZ(rho);

%compare within craniotomy to other areas
rho_in = cellfun(@(x) rho(x==1), within_mask,'UniformOutput',0);
rho_in = [rho_in{:}]; %merge sites
rho_in = rho_in(:);
%NaN out craniotomy sites
rho_out = rho;
for i = 1:numel(within_mask)
   rho_out(within_mask{i}==1)=NaN;
end
%remove NaN from both
rho_out(isnan(rho_out))=[];
rho_in(isnan(rho_in))=[];

end

%%
