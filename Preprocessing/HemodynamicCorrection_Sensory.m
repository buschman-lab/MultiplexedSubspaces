function [dff, dff_b, dff_v] = HemodynamicCorrection_Sensory(stack, opts)
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
    
    %Trim until both the same length (may be one frame discrepancy)
    min_length = cellfun(@(x) size(x,3), {stack_b,stack_v},'UniformOutput',0);
    stack_b = stack_b(:,:,1:min(cell2mat(min_length)));
    stack_v = stack_v(:,:,1:min(cell2mat(min_length)));
    
    %remove the masked pixels by setting to NaN
    masked_pxls = (stack_b==0);
    stack_b(masked_pxls)=NaN;
    stack_v(masked_pxls)=NaN;
   
    %pixelwise hemodynamic correction   
    [nX,nY,nZ] = size(stack_b);
    [stack_b, bad_col] = conditionDffMat(stack_b);
    [stack_v, ~] = conditionDffMat(stack_v);
    stack_v_corrected = NaN(size(stack_v));
    for i = 1:size(stack_v_corrected,2)  
        temp = stack_v(:,i); 
        %smooth with ~450ms gaussian
%         temp = smoothdata(temp,'gaussian',floor(opts.fps/2)); %opts.fps is per wavelength (if multiplexed) 
        %least square linear regression to find scaling factor and intercept        
        coef = polyfit(temp,stack_b(:,i),1); 
        correction = @(x) coef(1)*x+coef(2);
        stack_v_corrected(:,i) = correction(temp);
    end

    %remake into full size
    stack_b = conditionDffMat(stack_b,bad_col,[],[nX,nY,nZ]);    
    stack_v_corrected = conditionDffMat(stack_v_corrected,bad_col,[],[nX,nY,nZ]);
    
    %calculate the dff for each (as in wekselblatt, 2017)
    [dff_v,~] = makeDFF(stack_v_corrected, opts);
    [dff_b,~] = makeDFF(stack_b, opts);
    
    %Subtraction correction DFF = dff_b-dff_v;
    dff = dff_b-dff_v;       
       
end

    
        
% %     Smooth V to remove high frequency non hemodynamic activity
% %     
% %     
% %     least square linear regression to find scaling factor and intercept
% %     coef = polyfit(stack_v(~isnan(stack_v)),stack_b(~isnan(stack_b)),1);           
% %         
% %     correction equation    
% %     correction = @(x) coef(1)*x+coef(2);
% %     stack_v_corrected = reshape(correction(stack_v(:)),size(stack_v));  
% %     
% %     confirm that you scaled the correct direction
% %     err = nanmean(stack_v_corrected(:)-stack_b(:));
% %     err2 = nanmean(stack_v(:)-stack_b(:));
% %     
% %     calculate the dff for each    
% %     [dff_v,~] = makeDFF(stack_v_corrected, opts);
% %     [dff_b,~] = makeDFF(stack_b, opts);
% %     
% %     Subtraction correction DFF = dff_b-dff_v;
% %     dff = dff_b-dff_v;
% %        
