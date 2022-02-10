function [cstat,idx] = compareVectors(x,y,type,self_flag)
%Camden - timeless


idx = combvec([1:size(x,2)],[1:size(y,2)]); %pairwise comparison between all col in x and y

if self_flag
    idx(:,diff(idx)==0)=[];
end

cstat = NaN(size(idx,1),1);
for i = 1:size(idx,2)
    switch type
        case 'angle'            
            cstat(i) = AngleBetweenWeights(x(:,idx(1,i)),y(:,idx(2,i)),'none');
        case 'corr'
            cstat(i) = fisherZ(corr(x(:,idx(1,i))-nanmean(x(:,idx(1,i))),y(:,idx(2,i))-nanmean(y(:,idx(2,i)))));        
        case 'sse'
            cstat(i) = sum((x(:,idx(1,i))-y(:,idx(2,i))).^2);
        case 'mse'
            cstat(i) = nanmean((x(:,idx(1,i))-y(:,idx(2,i))).^2);            
        otherwise
            error('unknown consistency type');
    end
end





end %function 
