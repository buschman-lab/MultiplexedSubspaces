function [x_thresh, x_bin, x_onset] = ThresholdMatrix(x,val,method)
%x is motifs x time matrix or vector

x_bin = NaN(size(x)); 
for i = 1:size(x,1)
   temp = x(i,:);   
   switch method
       case 'std'
           x_bin(i,:) = temp>=(mean(temp)+(val * std(temp))); 
       case 'value'
           x_bin(i,:) = temp>=val;
       otherwise
           error('Unregonized method')
   end
end

%get the thresholded values
x_thresh = zeros(size(x));
x_thresh(x_bin==1) = x(x_bin==1); 

%1 = onset, -1 = offset. add back the first column  
x_onset = cat(2,zeros(size(x_bin,1),1),diff(x_bin==1,1,2));
x_onset = cellfun(@(x) find(x==1), num2cell(x_onset,2),'UniformOutput',0);
