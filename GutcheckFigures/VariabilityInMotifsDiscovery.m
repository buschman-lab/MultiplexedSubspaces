function VariabilityInMotifsDiscovery()

%% Fitting section AIC
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\Mouse9025-4-20-2018_DFFCombined_chunk2.mat');
opts = general_params_vpa;
opts.K
opts.repeat_fits = 5; 
num_runs = 5;
lambda = 0.003;
opts.fit_criterion = 'AIC';
%for reproducibility
rng('default');

%Fit Motifs To Training Data And Collect Statistics
idx = ones(num_runs,1);
crit = ones(num_runs,opts.repeat_fits);
W_temp = cell(num_runs,opts.repeat_fits);
H_temp = cell(num_runs,opts.repeat_fits);
stats_train_temp =cell(num_runs,opts.repeat_fits);
for cur_run = 1:num_runs
    for cur_fit = 1:opts.repeat_fits %fit multiple times due to random initialization
        if opts.verbose; fprintf('\nFitting Training Data Run %d of %d Fit %d of %d',cur_run,num_runs,cur_fit,opts.repeat_fits); end
        [W_temp{cur_run,cur_fit},H_temp{cur_run,cur_fit},stats_train_temp{cur_run,cur_fit}] = fpCNMF(data_train,'L',opts.L,'K',opts.K,'non_penalized_iter',...
            opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
            'speed','fast','verbose',opts.verbose,'lambda',lambda,...
            'ortho_H',opts.ortho_H,'w_update_iter',opts.w_update_iter,...
            'sparse_H',opts.sparse_H);  
        %Remove Empty Motifs 
        [W_temp{cur_run,cur_fit},H_temp{cur_run,cur_fit}] = RemoveEmptyMotifs(W_temp{cur_run,cur_fit},H_temp{cur_run,cur_fit});
    end


    %choose best fit
    [idx(cur_run), crit(cur_run,:)] = InternallyValidateWs(data_train,W_temp(cur_run,:),H_temp(cur_run,:),opts.fit_criterion,0);

end

%% Load the data 
load('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\GutcheckFigures\GutcheckData\VariabilityInMotifsDiscovery_firstfile_deconvolve_1sec_Mouse9025-4-20-2018_Chunk2.mat');
%% rerun for BIC
crit2 = ones(num_runs,opts.repeat_fits);
idx2 =ones(num_runs,1);
for cur_run = 1:num_runs    
    [idx2(cur_run), crit2(cur_run,:)] = InternallyValidateWs(data_train,W_temp(cur_run,:),H_temp(cur_run,:),'BIC',0);
end

%% Compute the dynamicness and spatiotemporal variabiility of the motifs 
dynamicness = NaN(num_runs,opts.repeat_fits);
stvar = NaN(num_runs,opts.repeat_fits);
for cur_run = 1:num_runs  
    for cur_fit = 1:opts.repeat_fits 
       temp = W_temp{cur_run,cur_fit};
       dynamicness(cur_run,cur_fit) = 1 - fisherInverse(nanmean(arrayfun(@(n) nanmean(fisherZ(corr(squeeze(temp(:,n,:)),nanmean(squeeze(temp(:,n,:)),2)))),1:size(temp,2),'UniformOutput',1)));   
       stvar(cur_run,cur_fit) = nanmean(arrayfun(@(n) nanvar(reshape(squeeze(temp(:,n,:)),1,numel(temp(:,n,:)))),1:size(temp,2),'UniformOutput',1));
    end
end


%% Compute the fit to the testing data
%Fit Motifs To Training Data And Collect Statistics
stats_test_temp =cell(num_runs,opts.repeat_fits);
for cur_run = 1:num_runs
    for cur_fit = 1:opts.repeat_fits %fit multiple times due to random initialization
        if opts.verbose; fprintf('\nFitting Test Data Run %d of %d Fit %d of %d',cur_run,num_runs,cur_fit,opts.repeat_fits); end       
        [~,~,stats_test_temp{cur_run,cur_fit}] = fpCNMF(data_test,'non_penalized_iter',...
            opts.non_penalized_iter,'penalized_iter',25,...
            'speed','fast','verbose',0,'lambda',lambda,'w_update_iter',0,...
            'ortho_H',opts.ortho_H,'W',W_temp{cur_run,cur_fit},'sparse_H',opts.sparse_H);
    end
end

%% Plotting Section 
col = {'k','b','k','b','k'};

%AIC Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    plot(x,crit(cur_run,:),'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),crit(cur_run,(idx(cur_run))),'marker','x','linewidth',2,'color','r','markersize',20)
end
title('Comparing AICr Across multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('AICr (AU)');

%BIC Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    plot(x,crit2(cur_run,:),'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx2(cur_run)),crit2(cur_run,(idx2(cur_run))),'marker','x','linewidth',2,'color','g','markersize',20)
end
title('Comparing BICr Across multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('BICr (AU)');

%PEV Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    pev = cellfun(@(x) x.pev, stats_train_temp(cur_run,:),'UniformOutput',1);
    plot(x,pev,'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),pev(idx(cur_run)),'marker','x','linewidth',2,'color','r','markersize',20)
    plot(x(idx2(cur_run)),pev(idx2(cur_run)),'marker','x','linewidth',2,'color','g','markersize',20)
end
title('Comparing PEV Across multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('PEV');

%PEV testing Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    pev = cellfun(@(x) x.pev, stats_test_temp(cur_run,:),'UniformOutput',1);
    plot(x,pev,'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),pev(idx(cur_run)),'marker','x','linewidth',2,'color','r','markersize',20)
    plot(x(idx2(cur_run)),pev(idx2(cur_run)),'marker','x','linewidth',2,'color','g','markersize',20)
end
title({'Comparing cross validation PEV','Across multiple Motif Fits/Runs'});
xlabel('Fits (colored by run');
ylabel('PEV');

% #Motifs Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    num = cellfun(@(x) x.n_motifs, stats_train_temp(cur_run,:),'UniformOutput',1);
    plot(x,num,'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),num(idx(cur_run)),'marker','x','linewidth',2,'color','r','markersize',20)
    plot(x(idx2(cur_run)),num(idx2(cur_run)),'marker','x','linewidth',2,'color','g','markersize',20)
end
title('Comparing #Motifs Across multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('#Motifs');

%Dynamicness Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    dyn = dynamicness(cur_run,:);
    plot(x,dyn,'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),dyn(idx(cur_run)),'marker','x','linewidth',2,'color','r','markersize',20)
    plot(x(idx2(cur_run)),dyn(idx2(cur_run)),'marker','x','linewidth',2,'color','g','markersize',20)
end
title('Dynamicsness Across Multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('average motif Dissimilarity(1-rho)');

%Spatiotemporal Variability Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    stv = stvar(cur_run,:);
    plot(x,stv,'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),stv(idx(cur_run)),'marker','x','linewidth',2,'color','r','markersize',20)
    plot(x(idx2(cur_run)),stv(idx2(cur_run)),'marker','x','linewidth',2,'color','g','markersize',20)
end
title({'Spatiotemporal Variability';'Across Multiple Motif Fits/Runs'});
xlabel('Fits (colored by run');
ylabel('average motif variability');

end %function
















