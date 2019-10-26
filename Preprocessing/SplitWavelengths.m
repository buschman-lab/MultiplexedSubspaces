function [stack_split] = SplitWavelengths(stack,pattern)

%repeat the wavelength strobe pattern to cover all frames
seq = repmat(pattern,1,ceil(size(stack,3)/length(pattern)));

%trim to exact length of recording
seq = seq(1:size(stack,3));

%split according to the unique wavelengths 
stack_split = arrayfun(@(x) stack(:,:,seq==x), unique(seq), 'UniformOutput',0);

end %function end





