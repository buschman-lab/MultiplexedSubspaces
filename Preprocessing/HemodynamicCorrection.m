function [dff, dff_b, dff_v] = HemodynamicCorrection(stack, opts)
    %Camden MacDowell 2019
    %Followed allen et al., 2017 neuron and Musall et al., 2019 Nature
    %Neuro subtraction method for hemodynamic correction. 

    fprintf('Performing hemodynamic correction');    
    stack = SplitWavelengths(stack,opts.wavelength_pattern);
    
    %Get hemo stack
    stack_v = stack{opts.correction_wavelength};  

    %neural signal
    stack(opts.correction_wavelength) = [];
    stack_b = stack{1};  
    
    %Trim until both the same length (may be one frame discrepancy
    min_length = cellfun(@(x) size(x,3), {stack_b,stack_v},'UniformOutput',0);
    stack_b = stack_b(:,:,1:min(cell2mat(min_length)));
    stack_v = stack_v(:,:,1:min(cell2mat(min_length)));
    
    %least square linear regression to find scaling factor and intercept
    coef = polyfit(stack_v(~isnan(stack_v)),stack_b(~isnan(stack_b)),1);
    
    %correction equation    
    correction = @(x) coef(1)*x+coef(2);
    stack_v_corrected = reshape(correction(stack_v(:)),size(stack_v));  
    
    %calculate the dff for each    
    [dff_v,~] = makeDFF(stack_v_corrected, opts);
    [dff_b,~] = makeDFF(stack_b, opts);
    
    %Subtraction correction DFF = dff_b-dff_v;
    dff = dff_b-dff_v;
       
end
