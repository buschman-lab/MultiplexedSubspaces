function [sig_motif,weight_motif,pref_motif,neu_area,discrim_mat, num_discrim] = SingleUnitANOVA(EphysPath,motif_fits,win)
%Camden MacDowell - timeless
%Compute the ANOVA over time for each neuron comparing activity to
%different motifs

if nargin <3; win = [-5 15]; end %post onset during with which to average over

%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe');
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%get the onsets for each motif
trig_st = cell(1,numel(motif_onset));
for i = 1:numel(motif_onset)
    [~,temp] = ParseByOnset([],st_norm,motif_onset,win,i);
    trig_st{i} = cat(1,temp{:});
end


%get the onsets for each motif
trig_st = cell(1,numel(motif_onset));
for i = 1:numel(motif_onset)
    [~,temp] = ParseByOnset([],st_norm,motif_onset,win,i);
    trig_st{i} = cat(1,temp{:});
end

%loop through each neuron
n = size(trig_st{1},1);
sig_motif = NaN(n,sum(abs(win))); %discriminatability (pval)
pref_motif = NaN(n,sum(abs(win))); %identity
weight_motif = NaN(n,sum(abs(win))); %weighting
unique_motif = cell(n,sum(abs(win))); %list of the pairs of motifs that can be disciminated between
for cur_n = 1:n
    %grab that neuron's activity in response to each motif
    neu_activity = cellfun(@(x) squeeze(x(cur_n,:,:)), trig_st,'UniformOutput',0);

    %loop through each point in time
    for cur_t = 1:size(neu_activity{1},1)
        cur_point = cellfun(@(x) x(cur_t,:),neu_activity,'UniformOutput',0);
        %assign motif 
        grp = arrayfun(@(n) n*ones(1,numel(cur_point{n})), 1:numel(cur_point),'UniformOutput',0);
        cur_point = cat(2,cur_point{:})';
        grp = cat(2,grp{:})';
        [sig_motif(cur_n,cur_t),~,stats] = anova1(cur_point,categorical(grp),'off'); %add categorical just to be sage
        %save off the 'preferential' motif
        [weight_motif(cur_n,cur_t),pref_motif(cur_n,cur_t)] = max(stats.means); %firing rate and identify of the preferred motif
        
        %also get the unique groups that are distinguishable
        [results,~] = multcompare(stats,'CType','tukey-kramer','display','off');        
        
        %get the pairs of motifs that are distinguishable by this neuron
        idx = find(results(:,end)<0.05);
        unique_motif{cur_n,cur_t} = results(idx,1:2); 
    end
end

%test similarity in motif preference across neurons
discrim_mat = NaN(n,n);
num_discrim = NaN(n,1);
for cur_t = 1:size(sig_motif,2)
    for i = 1:n
        for j = 1:n
           %how many of neurons motifs are also discriminatble in the target neuron
           discrim_mat(i,j,cur_t) = sum(arrayfun(@(n) sum(ismember(unique_motif{i,cur_t}(n,:),unique_motif{j,cur_t},'rows')), 1:size(unique_motif{i,cur_t},1),'UniformOutput',1)) / size(unique_motif{i,cur_t},1);       
           %NaN for any that were unable to discriminate anything
        end
        num_discrim(i,cur_t) = size(unique_motif{i,cur_t},1);
    end
end



% if verbose == 1       
%     figure('position',[684    68   339   928]); hold on;
%     %anything below our significance threshold (per timepoint) denote as NaN
%     sig_map = sig_motif; 
%     sig_map(sig_map>=0.05/numel(sig_map))=NaN;
%     sig_map(sig_map<0.05/numel(sig_map))=1;
%     
%     %colormap per motif
%     col = getColorPalet(15);
%     col = cat(1,[1 1 1],col);
%     temp = pref_motif;
%     temp(isnan(sig_map))=0;
%     imagesc(temp); colormap(col); hold on; 
%     ylim([0.5,size(sig_map,1)]);
%     set(gca,'ydir','reverse')
%     c = getColorPalet(30);
%     PlotProbeAnatomy(gca, neu_area, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:));
%     set(gca,'ytick','')
%     line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset   
%     
%     %loop through each motif
%     for i = 1:15
%         figure('position',[684    68   339   928]); hold on; 
%         temp = pref_motif;
%         temp(isnan(sig_map))=0; 
%         temp(temp~=i)=0;
%         imagesc(temp); colormap(cat(1,col(1,:),col(i,:))); hold on; 
%         ylim([0.5,size(sig_map,1)]);
%         set(gca,'ydir','reverse')
%         c = getColorPalet(30);
%         PlotProbeAnatomy(gca, neu_area, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:));
%         set(gca,'ytick','')
%         line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset  
%         title(sprintf('preference for motif %d',i),'fontweight','normal')
%     end
%     [~,rec_name] = fileparts(motif_fits);    
%     saveCurFigs(get(groot, 'Children'),{'-dpng'},[rec_name{1}(1:30)],save_dir,0); close all    
% end



end %function end























