function dff = HemodynamicCorrection(stack, opts)
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
    
    stack_h = NaN(size(temp));
    for pxl = 1:(nX*nY)
        stack_h(pxl,:) = filter(kern,1,temp(pxl,:));
    end
    
    %neural signal
    stack(opts.correction_wavelength) = [];
    stack_n = reshape(stack{1},[nX*nY,size(stack{1},3)]);    
    
    %Trim until both the same length (may be one frame discrepancy
    min_length = cellfun(@(x) size(x,2), {stack_n,stack_h},'UniformOutput',0);
    stack_n = stack_n(:,1:min(cell2mat(min_length)));
    stack_h = stack_h(:,1:min(cell2mat(min_length)));
    
    %least square linear regression to find scaling factor and intercept
    coef = polyfit(stack_h(:),stack_n(:),1);
    
    %correction equation    
    correction = @(x) coef(1)*x+coef(2);
    stack_h_corrected = reshape(correction(stack_h(:)),[nX*nY,size(stack_h,2)]);
    
    %reshape stacks to 3D (do this in case you ever use non-square stacks
    stack_h_corrected = reshape(stack_h_corrected,nX,nY,nZ);
    stack_n = reshape(stack_n,nX,nY,nZ);
    
    %calculate the dff for each    
    [dff_h,~] = makeDFF(stack_h_corrected, opts);
    [dff_n,~] = makeDFF(stack_n, opts);
    
    %Subtraction correction DFF = DFF_n-DFF_h;
    dff = dff_n-dff_h;
       
end
