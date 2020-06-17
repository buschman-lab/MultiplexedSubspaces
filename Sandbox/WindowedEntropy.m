function WindowedEntropy(H,dur)
%Windowed Motif Entrop

if nargin <2; dur =120*15; end 
H = H_orig;
%load the motifs and convolve H
motifs = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\AverageDPs_1.mat','W_clust_smooth');
motifs = motifs.W_clust_smooth;
M = size(H,1); %number of motifs

%reconvolve
for K = 1:size(H,1)
%     H(K,:) = tensor_convolve(nanmean(motifs(:,K,:),1),H(K,:));    
    H(K,:) = tensor_convolve(nansum(motifs(:,K,:),1),H(K,:));    %maybe try a sum... this would help deal with size differences
    %rolling integral to smooth it
    H(K,:) = movsum(H(K,:),15*5);
end

%smooth as well

%plot the variance over time
figure; hold on; 
H_var = movvar(H(:,60*15:end)',dur); %clip the first minute due to transiet at the begining
plot(H_var,'linewidth',1,'color',[0.75 0.75 0.75 1])
plot(nanmean(H_var,2),'linewidth',1,'color',[0 0 0 1])

%threshold the H values with moving standard deviation
H_thresh=H;
for i =1:size(H,1)   
   temp = H(i,:);
%    threshval = 2*movstd(temp,15*30,'omitnan');
   threshval = 2*std(temp,[],2);
   temp(temp<threshval)=0;
   H_thresh(i,:)=temp;  
end

%get the maximum motif at each timepoint 
[~,idx] = max(H_thresh,[],1);
bad_idx = sum(H_thresh,1); %remove indices where no motifs are activity (otherwise that will be identified as motif 1 activity)
idx(bad_idx==0)=NaN;

%window this index
dur = 15*15;
windowed_data=WindowData(idx,dur);

%for each window, compute the entropy of each motif
ent = NaN(numel(windowed_data),1);
all_pi =NaN(numel(windowed_data),M);
for w = 1:numel(windowed_data)
    %get probability of each motif 
    p_i = NaN(1,M+1);
    for m = 1:M %loop through motifs
       p_i(m) = sum(windowed_data{w}==m)/numel(windowed_data{w});         
       all_pi(w,m) = p_i(m);
    end
    p_i(M+1) = 1-nansum(p_i); %the probability that it's non of the motifs
    ent(w) = -nansum(p_i.*log(p_i));
end    

tot_p_i = NaN(1,M+1);
for m = 1:M %loop through motifs
   tot_p_i(m) = sum(idx==m)/numel(idx);            
end
tot_p_i(M+1) = 1-nansum(tot_p_i);
tot_ent = -nansum(tot_p_i.*log(tot_p_i));

rel_ent = ent/tot_ent;
range(rel_ent)

fp = fig_params;
%%
%plot the relative entropy
kernel = 5; 
figure; hold on; 
% plot(rel_ent,'linewidth',1.25,'color',[0.25 0.25 0.25 0.25]);
plot(movmean(rel_ent,kernel),'linewidth',2,'color',[0.25 0.25 0.25]);
ylim([0.9,1.1])
xlim([1,numel(windowed_data)]);
ylabel('Relative Entropy')
xlabel('Time (mins)');
set(gca,'xtick',linspace(0,numel(windowed_data),5),'xticklabel',linspace(0,60,5));
fp.FormatAxes(gca);
set(gca,'position',[0.2032    0.1696    0.7018    0.7554]);
grid on; 

%plot the entropy
figure; hold on; 
plot(ent,'linewidth',1.25,'color',[0.25 0.25 0.25 0.25]);
plot(movmean(ent,kernel),'linewidth',2,'color',[0.25 0.25 0.25]);
% fp.SetTitle(gca,'Entropy of Motif Activity Changes Over Time')
ylabel('Entropy (bits)')
xlabel('Time (minute)');
set(gca,'xtick',linspace(0,numel(windowed_data),5),'xticklabel',linspace(0,60,5));
fp.FormatAxes(gca);
set(gca,'position',[0.2032    0.1696    0.7018    0.7554]);
grid on; 
%%

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','entropy','Z:\Users\Camden\GrantsWithTim\Ro1_2020',1);
end


