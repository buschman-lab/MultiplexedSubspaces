function PlotLambdaSelection(file_list,parameter_class)
%Camden MacDowell - timeless
%loads the cost and regularization and averages across chunks

if nargin <2
   [file_list,~] = GrabFiles('\w*fit\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\MotifDiscovery'}); %select the preprocessed data (not the '_processed');
end

if nargin <3; parameter_class = 'general_params_corticaldynamics'; end 
    
gp = loadobj(feval(parameter_class)); 

data = cellfun(@(x) load(x,'stats_train'),file_list,'UniformOutput',0);
cost_norm = cellfun(@(x) x.stats_train.lambdafit_cost,data,'UniformOutput',0);
reg_norm = cellfun(@(x) x.stats_train.lambdafit_reg,data,'UniformOutput',0);
lambda = cellfun(@(x) x.stats_train.lambda,data,'UniformOutput',0);

cost_norm = cat(1,cost_norm{:});
reg_norm = cat(1,reg_norm{:});
lambda = cat(1,lambda{:});

fp = fig_params;
figure; hold on; 
rectangle('position',[nanmean(lambda)-std(lambda),0,2*std(lambda),1],'FaceColor',[0.25 0.25 0.25 0.25])
shadedErrorBar(gp.lambda,nanmean(cost_norm),std(cost_norm),'lineprops',{'color','r'})
shadedErrorBar(gp.lambda,nanmean(reg_norm),std(reg_norm),'lineprops',{'color','b'})
plot(gp.lambda,nanmean(cost_norm),'ro','linewidth',1.25); 
plot(gp.lambda,nanmean(reg_norm),'bo','linewidth',1.25); 
ylim([0 1])
set(gca,'xscale','log');
xlabel('Lambda'); ylabel('AU (std)');
title('Automated Lambda Selection','Fontweight',fp.font_weight,'Fontsize',fp.font_size)

fp = fig_params;
figure; hold on; 
rectangle('position',[nanmean(lambda)-sem(lambda),0,2*sem(lambda),1],'FaceColor',[0.25 0.25 0.25 0.25])
shadedErrorBar(gp.lambda,nanmean(cost_norm),sem(cost_norm),'lineprops',{'color','r'})
shadedErrorBar(gp.lambda,nanmean(reg_norm),sem(reg_norm),'lineprops',{'color','b'})
plot(gp.lambda,nanmean(cost_norm),'ro','linewidth',1.25); 
plot(gp.lambda,nanmean(reg_norm),'bo','linewidth',1.25); 
ylim([0 1])
set(gca,'xscale','log');
xlabel('Lambda'); ylabel('AU (sem)');
title('Automated Lambda Selection','Fontweight',fp.font_weight,'Fontsize',fp.font_size)

end