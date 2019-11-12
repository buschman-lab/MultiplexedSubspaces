function y = MakeCellsEqual(x,dim,pad_flag)
%makes all cells in x the same size along dim by pad or trim
desired_size = (cellfun(@(x) size(x,dim), x,'UniformOutput',0));

if pad_flag
    desired_size = max([desired_size{:}]);
    y= cellfun(@(v) [v, nan(1, desired_size-numel(v))], x, 'UniformOutput', false);
else
    desired_size = min([desired_size{:}]);
    if dim ==1
        y= cellfun(@(v) v(1:desired_size,:), x, 'UniformOutput', false);    
    elseif dim ==2
        y= cellfun(@(v) v(:,1:desired_size), x, 'UniformOutput', false);
    elseif dim ==3
        y= cellfun(@(v) v(:,:,1:desired_size), x, 'UniformOutput', false);
    elseif dim ==4
        y= cellfun(@(v) v(:,:,:,1:desired_size), x, 'UniformOutput', false);
    end
end
end %function end