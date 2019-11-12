function data_smooth = SpatialGaussian(data)

%Applies spatial gaussian filter to all data

%get nan values; 
idx = isnan(data);

%if nan values have already been removed, get the pixels with no varainces
%(e.g. nan)
if sum(idx(:))==0
    temp = reshape(data,[size(data,1)*size(data,2),size(data,3)]);
    nan_rows = nanvar(temp,[],2)<=eps;
end

%remove nanvalues for smoothing
data(isnan(data))=0; %remove nan so don't spread while smoothing

data_smooth = NaN(size(data));
for i = 1:size(data,3)
    data_smooth(:,:,i) = imgaussfilt(data(:,:,i),[1 1],'filterdomain','spatial');
end

data_smooth(idx)=0; %set the NaN's back to zero since they now have non-zero values

if sum(idx(:))==0
    temp = reshape(data,[size(data,1)*size(data,2),size(data,3)]);
    temp(nan_rows,:) = 0;
    data_smooth = reshape(temp,[size(data,1),size(data,2),size(data,3)]);
end

end
    