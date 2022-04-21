function CreateExampleSubspaceFigures(data_rrr,data_cca,savedir, rec_name)
%run through all pairs and create the summary figures (with labels) ... and
%the comparison with CCA figures. Also include in title, why a given
%subspace was not accepted. 
%data_rrr and data_cca are for One motif


% load example data
% data_rrr = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA\Mouse 331 Recording 1RRR_muaflag1_motif3.mat');
% data_cca = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA\Mouse 331 Recording 1CCA_muaflag1_motif3.mat','r','paired_areas');
% savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SubspacesTemp';

if ~exist(savedir); mkdir(savedir); end

n = size(data_rrr(1).paired_areas,1);

%get reasons for invalidity
[~,stats,bad_idx_all, label] = isSubspace(data_rrr, 0);

close all; 
motif = data_rrr(1).motif;
for i = 1:n        
    bad_idx = bad_idx_all{1};
    %plot the rrr summary figure; 
    plot_rrrSummary(stats.full_score{i},data_rrr(1).cvl_rrr{i},data_rrr(1).cvl_fa{i},ones(size(bad_idx,1),1)); 
    validstr = cellfun(@(x) [x ' '], label(bad_idx(i,:)),'UniformOutput',0);
    validstr = [validstr{:}];
    if isempty(validstr)
        title(gca,sprintf('M%d | Valid', motif));
    else
        title(gca,sprintf('M%d | invalid: %s', motif, validstr),'Interpreter','none');
    end
    ylabel(sprintf('Performance: %s to %s',data_rrr(1).area_label{data_rrr(1).paired_areas(i,:)}))  

    %plot the distribution of beta weights and cutoffs
    figure; hold on; 
    y = abs(data_rrr(1).rrr_B{i}(:,1)); 
%     y = abs(data_rrr(1).rrr_B{i}(:,1:stats.ss_dim(i))); 
%     c = getColorPalet(stats.ss_dim(i));
%     arrayfun(@(n) histogram(y(:,n),'facecolor',c(n,:),'facealpha',0.25,'edgealpha',0), 1:size(y,2),'UniformOutput',0)
    histogram(y,'facecolor',[0.25 0.25 0.25],'facealpha',0.5,'edgealpha',0)
    yvals = get(gca,'ylim');
    plot([prctile(y,75)+5*iqr(y),prctile(y,75)+5*iqr(y)],yvals,'linestyle','--','color','r');   
    ylabel('# Neurons'); xlabel('RRR beta weights')
    legend('data','5x IQR')        
    title(sprintf('%s to %s',data_rrr(1).area_label{data_rrr(1).paired_areas(i,:)}))  
end

%summary of why excluded
%remove self
bad_idx_all = bad_idx_all{1};
bad_idx_all = bad_idx_all(bad_idx_all(:,1)==0,2:end);
temp = sum(bad_idx_all)
label(2:end)


[ccaval,rrrval,rho,deltaD] = rrrCCAoverlap(data_rrr,data_cca,1); 

saveCurFigs(get(groot, 'Children'),{'-dpng'},['ExampleFigures',sprintf(' %s motif %d',rec_name, motif)],savedir,0); close all

end %function end

