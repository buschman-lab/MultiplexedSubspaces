function [delta_x, delta_y] = leaveOneOutCV(x,y,num_cvs)
%Camden MacDowell - timeless
%compute the impact of indviidual neurons on identified CVs. 
%inputs: 
%x and y are two tensors neurons x timepoits x trials. \
%num_cvs is the significant cvs as discovered by significantCVs

%concatentate across trials
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

%mean substracted
x = x-nanmean(x);
y = y-nanmean(y);

%run full cca
[~,~,r] = canoncorr(x,y); 
r = r.^2; 

%repeat after dropping each neuron out in both brain areas
idx = 1:size(x,2);
delta_x = NaN(1,size(x,2));
for i = 1:size(x,2)    
    [~,~,r_temp] = canoncorr(x(:,~ismember(idx,i)),y); 
    r_temp =r_temp.^2;
    delta_x(i) = sum(r(1:num_cvs)-r_temp(1:num_cvs));
end

idx = 1:size(y,2);
delta_y = NaN(1,size(y,2));
for i = 1:size(y,2)    
    [~,~,r_temp] = canoncorr(x,y(:,~ismember(idx,i))); 
    r_temp =r_temp.^2;
    delta_y(i) = sum(r(1:num_cvs)-r_temp(1:num_cvs));
end



end %function end















