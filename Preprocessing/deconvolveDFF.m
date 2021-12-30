function [data,nanpxs,x,y,z] = deconvolveDFF(data,gp)
%from Process and Split Data

%condition data and remove nan pxls
[x,y,z] = size(data);
[data,nanpxs] = conditionDffMat(data); %nanpxs are the same on each iteration so fine to overwrite

%deconvolve data
switch gp.w_deconvolution
    case 'filter_thresh'
        fprintf('\n\tFiltering and thresholding data')
        data = filterstack(data, opts.fps, gp.w_filter_freq, gp.w_filter_type, 1, 0);
        %Remove negative component of signal and find bursts as in MacDowell 2020. Legacy. %Deconvolution with non-negative output is preferred (maintains more data, comes with own assumptions). 
        for px = 1:size(data,2)
           temp = data(:,px);
           temp(temp<=nanmean(temp(:))+gp.w_nstd*std(temp(:))) = 0;
           data(:,px) = temp;
        end
        
        data(data<(nanmean(data(:))+gp.w_nstd*nanstd(data(:))))=0;
    case 'lucric'
        fprintf('\n\tPerforming a Lucy-Goosey Deconvolution (Lucy-Richardson)\n')
        for px = 1:size(data,2)
           data(:,px) = lucric(data(:,px),gp.d_gamma,gp.d_smooth,gp.d_kernel);
        end
    case 'only_filter'
        fprintf('\n\tFiltering data')
        data = filterstack(data, opts.fps, gp.w_filter_freq, gp.w_filter_type, 1, 0);
    case 'fNN_permouse' %feedforward neural network trained per animal
        fprintf('loading pretrained neural network');
        %load train nn(see GeneratefNN_spock.sh)        
        if ispc 
            load(['Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging' filesep gp.fit_nn_fn],'trained_opts');
            temp = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\restingstate_processed_fn.mat','dff_list');
            net_indx = find(strcmp(temp.dff_list,fn)==1);
        else
            load(['/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging' filesep gp.fit_nn_fn],'trained_opts');
            temp = load('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging/restingstate_processed_fn.mat','dff_list_bucket');
            net_indx = find(strcmp(temp.dff_list_bucket,fn)==1);            
        end       
        params = trained_opts{net_indx}.feedforwardparams;
        net =trained_opts{net_indx}.shallowfeedforward;
        %loop through each pixel
        for px = 1:size(data,2)   
            if mod(px,round(0.1*size(data,2) )) ==0
                fprintf('\t%g%% Complete\n', round(px./size(data,2)*100,2));
            end              
            xtemp = createRollingWindow(data(:,px), params.win)'; %t-n:t-1        
            stPred = net(xtemp)';                    
            stPred(stPred<0)=0; %convert purelin to poslin/relu
            %NaN pad to match timepoints
            data(:,px)=NaN;
            data(ceil(params.win/2):end-floor(params.win/2),px)=stPred;
        end
        %remove pad
        padidx = sum(isnan(data),2)==size(data,2); 
        data(padidx,:)=[];
        z =  size(data,1);
    case 'fNN_all' %feedforward neural network trained across mice
        fprintf('loading pretrained neural network');
        %load train nn(see GeneratefNN_spock.sh)        
        if ispc 
            load(['Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging' filesep gp.fit_nn_fn],'trained_opts');
        else
            load(['/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging' filesep gp.fit_nn_fn],'trained_opts');
        end       
        params = trained_opts{1}.feedforwardparams;
        net =trained_opts{1}.shallowfeedforward;
        %loop through each pixel
        for px = 1:size(data,2)   
            if mod(px,round(0.1*size(data,2) )) ==0
                fprintf('\t%g%% Complete\n', round(px./size(data,2)*100,2));
            end            
            xtemp = createRollingWindow(data(:,px), params.win)'; %t-n:t-1        
            stPred = net(xtemp)';        
             stPred(stPred<0)=0; %convert purelin to poslin/relu
            %NaN pad to match timepoints
            data(:,px)=NaN;
            data(ceil(params.win/2):end-floor(params.win/2),px)=stPred;
        end
        %remove pad
        padidx = sum(isnan(data),2)==size(data,2); 
        data(padidx,:)=[];
        z =  size(data,1);        
    case 'none'
        fprintf('doing nothing');
            
    otherwise
        error('unknown w_deconvolution option. Check general parameters'); 
end
