function [x, all_areas] = LoadVariable(data,type,area_name,tempval,normtype)
%Camden - timeless
%combines data across recordings and motifs to return a single variable

if nargin <4; tempval = []; end %just a burner variable
if nargin <5; normtype = 'mean'; end
all_areas = data{1}(1).area_label;


x = NaN(numel(data),size(data{1},2),numel(all_areas),2000);
for cur_rec = 1:numel(data) %loop through each recording
    for cur_motif = 1:size(data{cur_rec},2)%loop throuch each motif
        for cur_area = 1:numel(all_areas)            
            idx = strcmp(data{cur_rec}(cur_motif).area_label,all_areas{cur_area});  
            if sum(idx) >0
                switch type
                    case 'ridge_performance'
                        x(cur_rec,cur_motif,cur_area) = cellfun(@(x) nanmean(bestLambda(x,0)), data{cur_rec}(cur_motif).cvl_ridge(idx),'UniformOutput',1);
                    case 'number_neu'
                        x(cur_rec,cur_motif,cur_area) = size(data{cur_rec}(cur_motif).area_val{idx==1},1);
                    case 'engagement_zscore'                        
                        temp = GetLocalTrialAverageStrength(data{cur_rec}(cur_motif).area_val,data{cur_rec}(cur_motif).area_label,'max',1,normtype);
                        x(cur_rec,cur_motif,cur_area) = temp(idx);
                    case 'engagement_max'                        
                        temp = GetLocalTrialAverageStrength(data{cur_rec}(cur_motif).area_val,data{cur_rec}(cur_motif).area_label,'max',2,normtype);
                        x(cur_rec,cur_motif,cur_area) = temp(idx);
                    case 'engagement_minmax'                        
                        temp = GetLocalTrialAverageStrength(data{cur_rec}(cur_motif).area_val,data{cur_rec}(cur_motif).area_label,'max',3,normtype);                        
                        x(cur_rec,cur_motif,cur_area) = temp(idx);
                    case 'engagement'
                        temp = GetLocalTrialAverageStrength(data{cur_rec}(cur_motif).area_val,data{cur_rec}(cur_motif).area_label,'max',0,normtype);
                        x(cur_rec,cur_motif,cur_area) = temp(idx);
                    case 'ttt_variability'
                        temp = TrialToTrialVariability(data{cur_rec}(cur_motif).area_val,data{cur_rec}(cur_motif).area_label,'sum',normtype);
                        x(cur_rec,cur_motif,cur_area) = temp(idx);  
                    case 'ttt_activity_mean'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %in the predictors
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);
                        x(cur_rec,cur_motif,cur_area,1:numel(fr)) = fr;                           
                    case 'trial_mean'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %in predictors
                        fr = fr(:,3:end,:);
                        fr = nanmean(nanmean(fr,3),2);
                        x(cur_rec,cur_motif,cur_area,1:numel(fr)) = fr;
                      case 'trial_baseline'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %in predictors
                        fr = fr(:,1:2,:);
                        fr = nanmean(nanmean(fr,3),2);
                        x(cur_rec,cur_motif,cur_area,1:numel(fr)) = fr; 
                    case 'ttt_activity_mean_all'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %%in the dependents
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);
                        x(cur_rec,cur_motif,cur_area,1:numel(fr)) = fr;  
                    case 'ttt_variability_other'
                        temp = TrialToTrialVariability(data{cur_rec}(cur_motif).area_val,data{cur_rec}(cur_motif).area_label,'sum',normtype);
                        temp = {cat(1,temp{:})};
                        x(cur_rec,cur_motif,cur_area) = temp;                    
                    case 'rel_performance'
                        temp = nanmean(1-data{cur_rec}(cur_motif).cvl_rrr{idx})/nanmean(bestLambda(data{cur_rec}(cur_motif).cvl_ridge{idx},0));                        
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp; 
                    case 'rrr_auc'
                        temp = nanmean(1-data{cur_rec}(cur_motif).cvl_rrr{idx})/nanmean(bestLambda(data{cur_rec}(cur_motif).cvl_ridge{idx},0)); 
                        if isempty(tempval)
                            temp = [0,temp,1];                            
                        else
                            temp = [0,temp(1:tempval),1];
                        end                    
                        temp = trapz(temp/numel(temp));
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = 1-temp;                         
                    case 'rrr_performance'
                        temp = nanmean(1-data{cur_rec}(cur_motif).cvl_rrr{idx});                        
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp;  
                    case 'ntrial'
                        x(cur_rec,cur_motif,cur_area) =  size(data{cur_rec}(cur_motif).area_val{idx},3);                               
                    case 'rrr_dim'
                        if isempty(tempval) 
                            tempval = 0.8; 
                        end
                        temp = nanmean(1-data{cur_rec}(cur_motif).cvl_rrr{idx})/nanmean(bestLambda(data{cur_rec}(cur_motif).cvl_ridge{idx},0)); 
                        temp = find(temp>=tempval,1,'first');
                        if isempty(temp); temp=NaN; end
                        x(cur_rec,cur_motif,cur_area) = temp;                         
                    case 'cca_dim'
                        x(cur_rec,cur_motif,cur_area) = numel(data{cur_rec}(cur_motif).r{idx}); 
                    case 'cca_strength'
                        x(cur_rec,cur_motif,cur_area) = sum(data{cur_rec}(cur_motif).r{idx});
                    case 'num_inactive'
                        temp = data{cur_rec}(cur_motif).inactive_idx;
                        temp(cellfun(@(x) numel(x), temp)==1)=[];
                        x(cur_rec,cur_motif,cur_area) = sum(temp{idx});
                    case 'rrr_beta'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        %flip so most neurons positive
%                         if sum(temp(temp>0))<sum(temp(temp<0)) %strongest is weighted positively
                        if sum(temp>0)<sum(temp<0)
                            temp = temp*-1;
                        end
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp; 
                    case 'rrr_V'
                        temp = data{cur_rec}(cur_motif).rrr_V{idx}(:,tempval)';
%                         if sum(temp(temp>0))<sum(temp(temp<0)) %strongest is weighted positively   
                        %flip so most neurons positive
                        if sum(temp>0)<sum(temp<0) 
                            temp = temp*-1;
                        end
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp;
                    case 'rrr_beta_noflip'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp; 
                    case 'rrr_V_noflip'
                        temp = data{cur_rec}(cur_motif).rrr_V{idx}(:,tempval)';
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp;                           
                    case 'numberPredNeu'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';                      
                        x(cur_rec,cur_motif,cur_area) = numel(temp);
                    case 'numberTargetNeu'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); 
                        x(cur_rec,cur_motif,cur_area) = size(fr,1);
                    case 'rrr_synapticweight'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);   
                        if sum(temp>0)/numel(temp)<0.5 %adjuts directionality so majority are positive
                            temp = -1*temp;
                        end   
                        temp = temp.*fr';          
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp;   

                    case 'rrr_synapticweight_rev'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);  
                        if sum(temp>0)/numel(temp)<0.5
                            temp = -1*temp;
                        end                  
                        temp = temp.*fr';          
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp;  

                    case 'dom_neurons_num'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);                    
                        temp = abs(temp.*fr');    
                        [temp,~] = sort(temp,'descend');%sort                        
                        temp = temp/sum(temp);
                        temp = find(cumsum(temp)>=0.5,1,'first');
                        x(cur_rec,cur_motif,cur_area) = temp;   

                    case 'dom_neurons_auc'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);                    
                        temp = abs(temp.*fr');    
                        [temp,~] = sort(temp,'descend');%sort                        
                        temp = cumsum(temp/sum(temp));  
                        auc = trapz([0 temp])/numel([0,temp]);
                        x(cur_rec,cur_motif,cur_area) = 1-auc;  

                    case 'dom_neurons_aucrev'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);                    
                        temp = abs(temp.*fr');    
                        [temp,~] = sort(temp,'descend');%sort                        
                        temp = cumsum(temp/sum(temp));  
                        auc = trapz([0 temp])/numel([0,temp]);
                        x(cur_rec,cur_motif,cur_area) = 1-auc; 

                    case 'dom_neurons_label'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);                    
                        temp = abs(temp.*fr');    
                        [~,tempidx] = sort(temp,'descend');%sort                        
                        x(cur_rec,cur_motif,cur_area,1:numel(tempidx)) = tempidx;  
                    
                    case 'rrr_synapticweight_sorted'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);                       
                        temp = (temp.*fr');                          
                        [unique_temp,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas);
                        temp = arrayfun(@(n) sort(temp(temp_reordered==n),'descend'),unique_temp,'UniformOutput',0);
                        temp = cat(2,temp{:});
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp;                          
                        
                    case 'rrr_beta_max'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)'; 
                        %since Vt determines the predictive sign, flip so that strongest is weighted to positive
                        if sum(temp>0)/numel(temp)<0.5
                            temp = -1*temp;
                        end                        
                        
                        x(cur_rec,cur_motif,cur_area,1:size(temp,2)) = temp/(max(abs(temp(:))));  
                    case 'mean_fr' %mean firing rate of each neuron during the motif (trial average)
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = nanmean(nanmean(fr,3),2);
                        
                        x(cur_rec,cur_motif,cur_area,1:size(fr',2)) = fr'; 
                    case 'mean_fr_local' %mean firing rate of each neuron during the motif (trial average)
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = nanmean(nanmean(fr,3),2);
                        
                        x(cur_rec,cur_motif,cur_area,1:size(fr',2)) = fr';        
                 
                    case 'trialvar_fr' %mean firing rate of each neuron during the motif (trial average)
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);
                       
                        x(cur_rec,cur_motif,cur_area,1:size(fr',2)) = fr';   
                    case 'svca' %see SVCA
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = fr-nanmean(fr,3);
                        if ~isempty(tempval)
                            ndim = SVCA(fr,100,tempval);
                        else
                            ndim = SVCA(fr);
                        end
                        
                        x(cur_rec,cur_motif,cur_area) = ndim; 
                        
                    case 'svca_across' %see SVCA
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = fr-nanmean(fr,3);
                        if ~isempty(tempval)
                            [~,temp] = SVCA(fr,100,tempval);
                        else
                            [~,temp] = SVCA(fr);
                        end
                        
                        x(cur_rec,cur_motif,cur_area,1:size(temp,1)) = temp;                         

                    case 'pca'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = fr-nanmean(fr,3);
                        fr = reshape(fr,[size(fr,1),size(fr,2)*size(fr,3)]);
                        [~,~,~,~,ndim] = pca(fr');
                        if ~isempty(tempval)
                            ndim = find(cumsum(ndim)>tempval,1,'first');
                        else
                            ndim = find(cumsum(ndim)>80,1,'first');
                        end                        
                        x(cur_rec,cur_motif,cur_area) = ndim; 

                    case 'pca_auc'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = fr-nanmean(fr,3);
                        fr = reshape(fr,[size(fr,1),size(fr,2)*size(fr,3)]);
                        [~,~,~,~,ndim] = pca(fr');
                        ndim = cat(1,0,cumsum(ndim)/max(cumsum(ndim)));
                        ndim = trapz(ndim/numel(ndim));    
                        x(cur_rec,cur_motif,cur_area) = ndim; 

                    case 'pca_auc_time'
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        y = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas 
                        pev = PCAperTime(fr,y,tempval,cur_rec,cur_motif,idx,data);
                        temp = arrayfun(@(n) find(cumsum(pev(n,:))>=tempval,1,'first'),1:size(pev,1),'UniformOutput',1);
                        x(cur_rec,cur_motif,cur_area,1:numel(temp)) = temp; 

                    case 'svca_trialaverage' %see SVCA
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        if ~isempty(tempval)
                            ndim = SVCA(fr,100,tempval);
                        else
                            ndim = SVCA(fr);
                        end
                        
                        x(cur_rec,cur_motif,cur_area) = ndim;     
                    case 'svca_trialvar_acrosstime' %see SVCA
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);   
                        %uncomment below for trial to trail weighted
                        fr = fr-nanmean(fr,3);
                        ndim = NaN(1,size(fr,2));
                        for q = 1:size(fr,2)                            
                            if ~isempty(tempval)
%                                 ndim(q) = SVCA(fr,100,tempval,q);
                                ndim(q) = SVCA(fr(:,q,:),100,tempval);
                            else
%                                 ndim(q) = SVCA(fr,100,0,q);
                                ndim(q) = SVCA(fr(:,q,:));
                            end
                        end

                        x(cur_rec,cur_motif,cur_area,1:numel(ndim)) = ndim; 
                        
                    case 'rsq_acrosstime' %compute the fit of the model to each timepoint
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        y = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        rsq = ErrPerTime(fr,y,tempval,cur_rec,cur_motif,idx,data);
                        x(cur_rec,cur_motif,cur_area,1:numel(rsq)) = rsq;  
   
                    case 'auc_acrosstime' %compute the fit of the model to each timepoint
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        y = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas       
                        auc = AUCPerTime(fr,y,tempval,cur_rec,cur_motif,idx,data);
                        x(cur_rec,cur_motif,cur_area,1:numel(auc)) = 1-auc;  

                    case 'auc_acrosstime_rev' %compute the fit of the model to each timepoint
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        y = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas       
                        auc = AUCPerTime(fr,y,tempval,cur_rec,cur_motif,idx,data);
                        x(cur_rec,cur_motif,cur_area,1:numel(auc)) = 1-auc;  
                        
                    case 'rsq_acrosstime_rev' %compute the fit of the model to each timepoint
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        y = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        rsq = ErrPerTime(fr,y,tempval,cur_rec,cur_motif,idx,data);       
                        x(cur_rec,cur_motif,cur_area,1:numel(rsq)) = rsq;                          
                        
                    case 'trialvar_fr_local' %mean firing rate of each neuron during the motif (trial average)
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);
                        x(cur_rec,cur_motif,cur_area) = sum(fr); 
                        
                    case 'trialavg_timecourse' %
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); 
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);                     
                        fr = squeeze(nanmean(nanmean(fr,3),1));                        
                        x(cur_rec,cur_motif,cur_area,1:size(fr,2)) = fr; 

                    case 'trialavg_timecourse_other' %
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); 
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);                     
                        fr = squeeze(nanmean(nanmean(fr,3),1));                        
                        x(cur_rec,cur_motif,cur_area,1:size(fr,2)) = fr; 

                    case 'trialvar_timecourse' %
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1});
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);  
                        fr = abs(fr-nanmean(fr,3));
                        fr = squeeze(nanmean(nanmean(fr,3),1));                        
                        x(cur_rec,cur_motif,cur_area,1:size(fr,2)) = fr; 

                    case 'trialvar_timecourse_rev' %
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0});
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);  
                        fr = abs(fr-nanmean(fr,3));
                        fr = squeeze(nanmean(nanmean(fr,3),1));                        
                        x(cur_rec,cur_motif,cur_area,1:size(fr,2)) = fr; 

                    case 'synapticweight_sum'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
                        %uncomment for trial avg weighted
%                         fr = nanmean(nanmean(fr,3),2); 
                        %uncomment for trial to trail weighted
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);
                        if sum(temp>0)/numel(temp)<0.5 %flip so that most beta weights positive (since will be multiplied by Vt
                            temp = -1*temp;
                        end   

                        temp = temp.*fr';
                        [unique_temp,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas);
                        temp = arrayfun(@(n) nansum(temp(temp_reordered==n)),unique_temp);
                        
                        x(cur_rec,cur_motif,cur_area,unique_temp) = temp';
                        x(cur_rec,cur_motif,cur_area,~ismember(1:numel(all_areas),unique_temp))=0; %remove the shared nan                      
                    case 'synapticweight_sum_abs'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
%                         fr = nanmean(nanmean(fr,3),2);
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);
                        temp = abs(temp.*fr'); 
                        [unique_temp,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas);
                        temp = arrayfun(@(n) nansum(temp(temp_reordered==n)),unique_temp);
                        
                        x(cur_rec,cur_motif,cur_area,unique_temp) = temp';
                        x(cur_rec,cur_motif,cur_area,~ismember(1:numel(all_areas),unique_temp))=0; %remove the shared nan
                    case 'synapticweight_mean_abs'
                        temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)';
                        fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
                        fr = normalizeToBaseline(fr,[1:2],normtype);
                        fr = fr(:,3:end,:);
%                         fr = nanmean(nanmean(fr,3),2);
                        fr = abs(fr-nanmean(fr,3));
                        fr = nanmean(nanmean(fr,3),2);
                        temp = abs(temp.*fr'); 
                        [unique_temp,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas);
                        temp = arrayfun(@(n) nanmean(temp(temp_reordered==n)),unique_temp);
                        
                        x(cur_rec,cur_motif,cur_area,unique_temp) = temp';
                        x(cur_rec,cur_motif,cur_area,~ismember(1:numel(all_areas),unique_temp))=0; %remove the shared nan
                    case 'rrr_bavg_max'
                        beta_temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)'; 
                        %since Vt determines the predictive sign, flip so that strongest is weighted to positive
                        if  sum(beta_temp>0)/numel(beta_temp)<0.5
                            beta_temp = -1*beta_temp;
                        end
                        beta_temp=beta_temp/(max(abs(beta_temp(:))));                        
                        beta_temp=beta_temp/(max(abs(beta_temp(:))));
                        [unique_temp,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas);
                        beta_temp = arrayfun(@(n) nanmean(beta_temp(temp_reordered==n)),unique_temp);
                        
                        x(cur_rec,cur_motif,cur_area,unique_temp) = beta_temp'; 
                        x(cur_rec,cur_motif,cur_area,~ismember(1:numel(all_areas),unique_temp))=0; %remove the shared nan
                    case 'rrr_bsum_max'
                        beta_temp = data{cur_rec}(cur_motif).rrr_B{idx}(:,tempval)'; 
                        %since Vt determines the predictive sign, flip so that strongest is weighted to positive
                        if sum(beta_temp>0)/numel(beta_temp)<0.5
                            beta_temp = -1*beta_temp;
                        end
                        beta_temp=beta_temp/(max(abs(beta_temp(:))));

                        [unique_temp,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas);                
                        beta_temp = arrayfun(@(n) nansum(beta_temp(temp_reordered==n)),unique_temp);                          
                        
                        x(cur_rec,cur_motif,cur_area,unique_temp) = beta_temp';  
                        x(cur_rec,cur_motif,cur_area,~ismember(1:numel(all_areas),unique_temp))=0; %remove the shared nan

                    case 'beta_region'                         
                        [~,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas);
                        x(cur_rec,cur_motif,cur_area,1:size(temp_reordered',2))=temp_reordered';
                    otherwise
                        error('unknown metric');
                end %switch
            else
                x(cur_rec,cur_motif,cur_area)=NaN;
            end %if/else
        end %area loop
    end %motif loop
end %rec loop

%remove extraneous 4th dimension
%get the one with the fewest nan
temp =reshape(x,[size(x,1)*size(x,2)*size(x,3),size(x,4)]);
x(:,:,:,(1+max(sum(~isnan(temp),2))):end)=[];

if ~isempty(area_name)
   x = squeeze(x(:,:,strcmp(all_areas,area_name),:));
end

end %function

function [unique_temp,temp_reordered] = reassign(data,cur_rec,cur_motif,idx,all_areas)
    %for recordings without all regions, need to reassign
    lbl = cellfun(@(x) find(strcmp(x,data{cur_rec}(cur_motif).area_label)), all_areas,'UniformOutput',0);
    temp = data{cur_rec}(cur_motif).grouping{idx};
    lbl(cellfun(@isempty,lbl))={NaN};
    lbl = [lbl{:}];
    temp_reordered = temp;
    for q = 1:numel(lbl)
        temp_reordered(temp==lbl(q)) = q;
    end                      
    %average beta per area
    unique_temp = unique(temp_reordered);
end %sub function 

function rsq = ErrPerTime(fr,y,tempval,cur_rec,cur_motif,idx,data)
    fr = normalizeToBaseline(fr,[1:2],normtype);
    y = normalizeToBaseline(y,[1:2],'mean');                        
    fr = fr(:,3:end,:);   
    y = y(:,3:end,:);
    fr = fr-nanmean(fr,3);
    y = y-nanmean(y,3);
    %load the betas and the V
    tempval = min(tempval,size(y,1));
    BB = data{cur_rec}(cur_motif).rrr_B{idx}(:,1:tempval);
    VV = data{cur_rec}(cur_motif).rrr_V{idx}(:,1:tempval);  
    rsq = NaN(1,size(fr,2));
    for q = 1:size(fr,2)      
       xx = fr(:,q,:);
       yy = y(:,q,:);
       xx = reshape(xx,[size(xx,1),size(xx,2)*size(xx,3)])';
       yy = reshape(yy,[size(yy,1),size(yy,2)*size(yy,3)])';
       Yhat = xx*BB*VV';
       rsq(q) = 1-NormalizedSquaredError(yy, Yhat);
    end 
end 

function auc = AUCPerTime(fr,y,tempval,cur_rec,cur_motif,idx,data)
    ndim = tempval;
    rsq = arrayfun(@(n) ErrPerTime(fr,y,n,cur_rec,cur_motif,idx,data), 1:ndim,'UniformOutput',0);
    rsq = cat(1,rsq{:});
    %normalize by the max so between 0 and 1
    rsq = rsq./max(rsq,[],1);
    auc = arrayfun(@(n) find(rsq(:,n)>0.8,1,'first'), 1:size(rsq,2));
    %get AUC
%     auc = arrayfun(@(n) trapz([0 rsq(:,n)'])/(1+ndim), 1:size(rsq,2));
end %funciotn 


function pev = PCAperTime(fr,y,tempval,cur_rec,cur_motif,idx,data)
    fr = normalizeToBaseline(fr,[1:2],normtype);
    fr = fr(:,3:end,:);   
    fr = fr-nanmean(fr,3);
	x = reshape(fr,[size(fr,1),size(fr,2)*size(fr,3)]);
    %loop through each PCA dimension and get the explained variance per
    %timepoint
    t = size(fr,2);
    n = size(fr,1)-1;
%     pev =[];
%     warning off
%     for q = 1:t  
%        yy = x(:,q:t:end);
%        [~,~,~,~,pev(q,:)] = pca(yy);
%     end


    %compute PCA on all data
    [coef, score, ~, ~, ~, mu] = pca(x);

    %loop through each PCA dimension and get the explained variance per
    %timepoint
    t = size(fr,2);
    n = size(score,2);
    pev = NaN(t,n);
    for q = 1:t  
        for n = 1:size(fr,1)-1
            Yhat = score(:,n)*coef(q:t:end,n)'+repmat(mu(q:t:end),size(score,1),1);
            yy = fr(:,q,:);
            yy = reshape(yy,[size(yy,1),size(yy,2)*size(yy,3)]); 
            pev(q,n) = PercentExplainedVariance(yy, Yhat, 0,0);
        end    
    end
end 


















