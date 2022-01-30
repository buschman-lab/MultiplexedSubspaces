function [dff_probe,offset] = LoadInsertionSiteDFF(ImgPath,ImgProbeLoc)

%load imaging 
load(ImgPath,'data_norm','nanpxs');  
data = conditionDffMat(data_norm',nanpxs); 
probe_loc = load(ImgProbeLoc); probe_loc = probe_loc.probe_coords;

x = sqrt(size(data_norm,1)+numel(nanpxs)); 
r = 2;
offset = {[2,0],[0,-1],[1,-1],[0 0]};
%parse the radius around the probes
dff_probe = NaN(size(data,3),numel(probe_loc));
for cur_probe = 1:numel(probe_loc)
    %create mask of radius around probe tip
    temp = probe_loc{cur_probe};
    mask = zeros(x,x);  
    temp(1,1)=temp(1,1)+offset{cur_probe}(1);
    temp(1,2)=temp(1,2)+offset{cur_probe}(2);
    mask(temp(1,2)-r:temp(1,2)+r,temp(1,1)-r:temp(1,1)+r)=1;
    %flatten
    mask = mask(:);
    mask(nanpxs)=[];        
    
    dff_probe(:,cur_probe) = nanmean(data_norm(mask==1,:),1); 
end    
