function windowed_data=WindowData(data,window)
%Camden MacDowell - timeless
%breaks data (N x T) into equal chunks of N x dur. 

[~,T] = size(data);

num_chunks = floor(T/window);

%remove the remainder and reshape into chunks
data_trim = data(:,1:end-mod(T,num_chunks*window));
data_trim = reshape(data_trim,[size(data_trim,1),num_chunks,window]);

windowed_data = arrayfun(@(x) squeeze(data_trim(:,x,:)), 1:num_chunks,'UniformOutput',0);

end %function end


