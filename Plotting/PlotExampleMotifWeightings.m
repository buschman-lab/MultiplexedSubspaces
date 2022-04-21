function PlotExampleMotifWeightings(cur_rec,motif_num,tp)
%camden macdowell
if nargin <3; tp = [100,100+(60*15*1)]; end

%reasonable timepoints = 

fp = fig_params_cortdynamics;

[~,~,~,~,motif_fits] = LoadDataDirectories(cur_rec);

%IMPORTANT: motif_fits returns motif fits in order test-->train. But true
%data is train-->test. CompileMotifOnsets fixes this, but any time that you
%load motif fits you'll want to reorder (like below). 
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%remove noise motif
motif_onset(fp.noisemotif) = [];

%%
%plot the motif weightings
col = repmat([0.75 0.75 0.75],numel(motif_onset),1); %grays
col(motif_num,:) = cat(1,fp.c_ff,fp.c_glm,fp.c_lr);

figure; hold on; 
for i = 1:size(col,1)
    if sum(ismember(motif_num,i))>0
        arrayfun(@(n) plot([motif_onset{i}(n),motif_onset{i}(n)],[0 1],'linestyle','-','color',col(i,:),'linewidth',2), 1:numel(motif_onset{i}));
    else
        arrayfun(@(n) plot([motif_onset{i}(n),motif_onset{i}(n)],[0 1],'linestyle','-','color',col(i,:),'linewidth',1), 1:numel(motif_onset{i}));
    end
end
xlim(tp)
% set(gca,'xtick',0:(15*10):(tp(2)-tp(1)),'xticklabel',(0:(15*10):(tp(2)-tp(1)))/15)
% xlabel('time (seconds)');
% ylabel('neurons')
fp.FormatAxes(gca); box on

% %plot the timecourse of the motif here 
% %only grab up to as many ass needed to get through tp
% data = load(motif_fits{i},'w','H','stats_refit');
% nn = ceil(tp(2)/size(data.H,2));
% weight_all = cell(1,nn);
% for i = 1:nn
%     data = load(motif_fits{i},'w','H','stats_refit');
%     weight = arrayfun(@(n) tensor_convolve(nanmax(data.w(:,n,:),[],1),data.H(n,:)),1:size(data.w,2),'UniformOutput',0);
%     weight_all{i} = cat(1,weight{:});
% %     weight_all{i} = data.stats_refit.rho_frame_per_motif;
% end
% weight = [weight_all{:}];
% weight(fp.noisemotif,:)=[];
% 
% %plot the motif weightings
% col = repmat([0.75 0.75 0.75],numel(motif_onset),1); %grays
% col(motif_num,:) = cat(1,fp.c_ff,fp.c_glm,fp.c_lr);


% figure; hold on; 
% for i = 1:size(col,1)
%     %get the one second after the motif is on
%     temp = weight(i,1:tp(2));
%     temp = temp-nanmean(temp);
%     temp(temp<0)=0;
% %     tempy = find(y(motif_num(i),:)==1);
%     tempy = find(y(i,:)==1);
%     if ~isempty(tempy)
%         for j = 1:numel(tempy)            
%             plot((tempy(j)):(tempy(j)+dur),temp((tempy(j)):(tempy(j)+dur)),'linewidth',2,'color',col(i,:))
%         end
%     end   
% end
%%


end %function end
