function area_val = RemoveEdgeTrials(area_val)
%Camden - timeless
%Removes the small edge case of trials that occur immediately at the start
%of a recording and so have no 'baseline' to substract. 

%all rows are going to be the same
atstart = sum(isnan(area_val{1}(1,:)),'all');

if atstart>0 %yes edge case | remove first trial
   area_val = cellfun(@(x) x(:,:,2:end), area_val,'UniformOutput',0); 
end


end