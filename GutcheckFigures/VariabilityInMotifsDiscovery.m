function VariabilityInMotifsDiscovery()

%% Fitting section
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\Mouse9025-4-20-2018_DFFCombined_chunk2.mat');
opts = general_params_vpa;
opts.K
opts.repeat_fits = 5; 
num_runs = 5;
lambda = 0.003;
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

%% Plotting Section 
col = {'k','b','k','b','k'};

%AIC Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    plot(x,crit(cur_run,:),'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),crit(cur_run,(idx(cur_run))),'marker','x','linewidth',2,'color','r')
end
title('Comparing AICr Across multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('AICr (AU)');

%PEV Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    pev = cellfun(@(x) x.pev, stats_train_temp(cur_run,:),'UniformOutput',1);
    plot(x,pev,'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),pev(idx(cur_run)),'marker','x','linewidth',2,'color','r')
end
title('Comparing PEV Across multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('PEV');

% #Motifs Figure
figure;  hold on 
for cur_run = 1:num_runs
    x = (5*(cur_run-1))+(1:opts.repeat_fits);
    num = cellfun(@(x) x.n_motifs, stats_train_temp(cur_run,:),'UniformOutput',1);
    plot(x,num,'marker','o','linewidth',2,'color',col{cur_run})
    plot(x(idx(cur_run)),num(idx(cur_run)),'marker','x','linewidth',2,'color','r')
end
title('Comparing #Motifs Across multiple Motif Fits/Runs');
xlabel('Fits (colored by run');
ylabel('#Motifs');

end %function
















