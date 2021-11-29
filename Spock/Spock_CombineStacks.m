function Spock_CombineStacks(folder_path,save_fn,parameter_class)
    if ~ispc
        addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    end
    
    warning('Camden you need to make sure spock is working for this code'); 
%      try       
    [stack, opts] = CombineStacks(folder_path,parameter_class);
    
    %legacy compatibility
    if ~isfield(opts,'savecompressed')
        opts.savecompressed=0;
    end
    
    %perform hemodynamic correction and make dff
    if numel(unique(opts.wavelength_pattern))>1 %if multiple wavelengths used        
       [dff, dff_b, ~] = HemodynamicCorrection(stack, opts); 
       %ImpactOfHemoCorrection(dff,dff_b,dff_h)
    else %Camden: do NOT add filterstack here... gums up the deconvolution with GLM (GLM kernel is then predominately the filter since the spiking is not filtered). 
       fprintf('\n No hemodynamic correction');             
       dff_b = makeDFF(stack, opts); 
       dff = [];
    end

    fprintf('\n Saving data');
    %Save off corrected data if available
    if isempty(dff)
       dff = dff_b;   
    end
    
    %flatten and compress
    if opts.savecompressed == 1 %full length, but flattened   
       [dff,nanpxs] = conditionDffMat(dff);
       save(save_fn,'dff','opts','nanpxs','-v7.3');        
    elseif opts.savecompressed == 2 %chunked, and flattened
       [dff,nanpxs] = conditionDffMat(dff); 
    else %full tensor
       save(save_fn,'dff','opts','-v7.3');
    end     

    %save off the uncorrected if desired
    if opts.save_uncorrected
        [path, fn] = fileparts(save_fn);
        dff = dff_b; %need to have it still named dff. 
        fn = [path filesep fn '_dff_uncorrected.mat'];
        %flatten and compress
        if opts.savecompressed    
           [dff,nanpxs] = conditionDffMat(dff);
           save(fn,'dff','opts','nanpxs','-v7.3');          
        else
           save(fn,'dff','opts','-v7.3');
        end          
    end
        
        %make figures
%         fh = DetectDictalEvents(save_fn);
%     catch
%         gp = general_params; 
%         warning('\nNo dffs with suffix %s were combined in %s',gp.stack_suffix, folder_path);
%     end
end

