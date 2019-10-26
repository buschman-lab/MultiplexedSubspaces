function [dff, dff_b, dff_v] = HemodynamicCorrection(stack, opts)
    %Camden MacDowell 2019
    %Followed allen et al., 2017 neuron and Musall et al., 2019 Nature
    %Neuro subtraction method for hemodynamic correction. 

    fprintf('Performing hemodynamic correction');    
    stack = SplitWavelengths(stack,opts.wavelength_pattern);

    %smooth the hemo stack per pixel 
    temp = stack{opts.correction_wavelength};   
    [nX,nY,nZ] = size(temp); 
    temp = reshape(temp,nX*nY,nZ);
    
    %400ms gaussian kernel
    kern = gausswin(round(400/(1000/opts.fps),0));
    
    stack_v = NaN(size(temp));
    for pxl = 1:(nX*nY)
        stack_v(pxl,:) = filter(kern,1,temp(pxl,:));
    end
    
    %neural signal
    stack(opts.correction_wavelength) = [];
    stack_b = reshape(stack{1},[nX*nY,size(stack{1},3)]);    
    
    %Trim until both the same length (may be one frame discrepancy
    min_length = cellfun(@(x) size(x,2), {stack_b,stack_v},'UniformOutput',0);
    stack_b = stack_b(:,1:min(cell2mat(min_length)));
    stack_v = stack_v(:,1:min(cell2mat(min_length)));
    
    %least square linear regression to find scaling factor and intercept
    coef = polyfit(stack_v(:),stack_b(:),1);
    
    %correction equation    
    correction = @(x) coef(1)*x+coef(2);
    stack_v_corrected = reshape(correction(stack_v(:)),[nX*nY,size(stack_v,2)]);
    
    %reshape stacks to 3D (do this in case you ever use non-square stacks
    stack_v_corrected = reshape(stack_v_corrected,nX,nY,nZ);
    stack_b = reshape(stack_b,nX,nY,nZ);
    
    %calculate the dff for each    
    [dff_v,~] = makeDFF(stack_v_corrected, opts);
    [dff_b,~] = makeDFF(stack_b, opts);
    
    %Subtraction correction DFF = dff_b-dff_v;
    dff = dff_b-dff_v;
       
end
