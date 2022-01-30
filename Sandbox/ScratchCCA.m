%scratch CCA
%%
% if cur_rec ==1
%     EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
%     [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% elseif cur_rec ==2
%     EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat';
%     [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_12_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% elseif cur_rec ==3
%     EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
%     [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% elseif cur_rec==4
%     EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_08_2021_RestingState_g0\ap_opts.mat';
%     [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_08_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% elseif cur_rec==5
%     EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
%     [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% elseif cur_rec==6
%     EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_10_2021_RestingState_g1\ap_opts.mat';
%     [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% end

%target matrix is a neruons x [trials] then you'll do xfold validation on


%Load the ephys for motif 9


win = [-5 15]; %post onset during with which to average over

EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
[motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});

%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe'); 
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%parse onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,win,7);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');


%% Question 1: Do significant subspaces exist? 
% Method: do CCA on two brain regions that the motif traverses. Significance is determined through trial permutation of activity in one of the brain regions. 
% Control: test whether cortical areas not implicated in the motif have significant subspaces. 
% Question 2: How much of the variance in local activity is captured by the communication subspace? 
% Method: Compute PCA on the activity of both areas in the subspace and normalize both the 
% CCA and PCA weightings across neurons to 0 and 1. Compute the angle between them and get the bootstrapped distribution. 
% Question 2: How much of the variance in local activity is captured by the communication subspace?
%motif 3; which includes both MOs and SSp-bfd 
compare_area = {'MOs','SSp-bfd'}; %VISp
dur = sum(abs(win))+1;
x = area_val{strcmp(area_label,compare_area{1}),:};
y = area_val{strcmp(area_label,compare_area{2}),:};
[a,b,U,V,r] = significantCVs(x,y,0.05,0);

if ~isempty(r); fprintf('\n\tsignificant subspace found between %s and %s for motif %d',compare_area{1},compare_area{2},3); end

%concatentate across trials
x_temp = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y_temp = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

%mean substracted
x_temp = x_temp-nanmean(x_temp);
y_temp = y_temp-nanmean(y_temp);
x_coef = pca(x_temp);
y_coef = pca(y_temp);

%get the angle between the subspace and pca weightings
cv_num = 1; 
cca_norm = a(:,cv_num)/max(a(:,cv_num));
pca_norm = x_coef(:,cv_num)/max(x_coef(:,cv_num)); 
x_theta = acosd( dot(pca_norm,cca_norm) / ( norm(pca_norm)*norm(cca_norm) ) );
%generate example figure
figure; hold on; 
bar(pca_norm,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor','r'); 
bar(cca_norm,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor','b'); 
cca_norm = b(:,cv_num)/max(b(:,cv_num));
pca_norm = y_coef(:,cv_num)/max(y_coef(:,cv_num)); 
y_theta = acosd( dot(pca_norm,cca_norm) / ( norm(pca_norm)*norm(cca_norm) ) );
%generate example figure
figure; hold on; 
bar(pca_norm,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor','r'); 
bar(cca_norm,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor','b'); 

%compute how much local variance is being captured by the subspace 
% [a_full,b_full] = canoncorr(x_temp,y_temp); 
% resid = x_temp - (x_temp*a_full(:,cv_num));
% exp_var = 1 - nanvar(resid(:))./nanvar(x_temp(:));
% resid = y_temp - (y_temp*b_full(:,cv_num));
% exp_var = 1 - nanvar(resid(:))./nanvar(y_temp(:));


%Question 3. Are regions sharing subspaces or are they different? 
%run CCA on all pairs of regions
area_label(strcmp(area_label,'cc'))=[];
area_label(strcmp(area_label,'na'))=[];
area_label(strcmp(area_label,'fiber tracts'))=[];

% compare_area = nchoosek(1:numel(area_label),2); 
% %ignore self and reverse comparisons
%let's start with two different pairs with MOs; BFD, MED, VISa. 
dur = sum(abs(win))+1;
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'MED'),:};
tic; [aa,bb,UU,VV,rr] = significantCVs(x,y,0.05,0); toc

x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'SSp-n'),:};
[aaa,bbb,UUU,VVV,rrr] = significantCVs(x,y,0.05,0);

x = area_val{strcmp(area_label,'SSp-bfd'),:};
y = area_val{strcmp(area_label,'SSp-n'),:};
[c,d,UUU,VVV,rrr] = significantCVs(x,y,0.05,0);
%so then compare bbb to d vs bbb' to d and b to c vs b' to c

%compare the angle between the same source area and each of its subspaces
cv_num=1;
cca_norm = a(:,cv_num)/max(a(:,cv_num));
temp = aa(:,cv_num)/max(aa(:,cv_num));
y_theta = acosd( dot(temp,cca_norm) / ( norm(temp)*norm(cca_norm) ) );

cca_norm = a(:,cv_num)/max(a(:,cv_num));
temp = aaa(:,cv_num)/max(aaa(:,cv_num));
y_theta = acosd( dot(temp,cca_norm) / ( norm(temp)*norm(cca_norm) ) );


%get the correlation in the beta weights of the subspaces
cv_num=1;
rho = corr(a(:,cv_num),aa(:,cv_num));
rho = corr(a(:,cv_num),aaa(:,cv_num));

%compute the permuted angle
cv_num=1;
cca_norm = a(:,cv_num)/max(a(:,cv_num));
temp = aaa(:,cv_num)/max(aaa(:,cv_num));
y_theta = acosd( dot(temp,cca_norm) / ( norm(temp)*norm(cca_norm) ) );
if y_theta>90
   y_theta = 180-y_theta; 
end

a_perm = arrayfun(@(n) a(randperm(size(a,1),size(a,1)),cv_num), 1:1000,'UniformOutput',0);
a_perm = cat(2,a_perm{:});
y_theta_perm = NaN(1,1000);
for i = 1:size(a_perm,2)
    cca_norm = a(:,cv_num)/max(a(:,cv_num));
    temp = a_perm(:,i)/max(a_perm(:,i));
    y_theta_perm(i) = acosd( dot(temp,cca_norm) / ( norm(temp)*norm(cca_norm) ) );
    if y_theta_perm(i)>90
       y_theta_perm(i) = 180-y_theta_perm(i); 
    end 
end

figure; hold on; histogram(y_theta_perm); plot([y_theta,y_theta],[0 100],'r');

%compute the permuted angle
cv_num=1;
cca_norm = a(:,cv_num)/max(a(:,cv_num));
temp = aa(:,cv_num)/max(aa(:,cv_num));
y_theta = acosd( dot(temp,cca_norm) / ( norm(temp)*norm(cca_norm) ) );
if y_theta>90
   y_theta = 180-y_theta; 
end

a_perm = arrayfun(@(n) a(randperm(size(a,1),size(a,1)),cv_num), 1:1000,'UniformOutput',0);
a_perm = cat(2,a_perm{:});
y_theta_perm = NaN(1,1000);
for i = 1:size(a_perm,2)
    cca_norm = a(:,cv_num)/max(a(:,cv_num));
    temp = a_perm(:,i)/max(a_perm(:,i));
    y_theta_perm(i) = acosd( dot(temp,cca_norm) / ( norm(temp)*norm(cca_norm) ) );
    if y_theta_perm(i)>90
       y_theta_perm(i) = 180-y_theta_perm(i); 
    end 
end
figure; hold on; histogram(y_theta_perm); plot([y_theta,y_theta],[0 100],'r');

%compute the permuted correlation
cv_num=1;
rho = abs(corr(a(:,cv_num),aaa(:,cv_num)));
rho_perm = arrayfun(@(n) abs(corr(a(:,cv_num),a_perm(:,n))),1:size(a_perm,2));

a_perm = arrayfun(@(n) a(randperm(size(a,1),size(a,1)),cv_num), 1:1000,'UniformOutput',0);
a_perm = cat(2,a_perm{:});
rho_perm2 = arrayfun(@(n) abs(corr(a(:,cv_num),a_perm(:,n))),1:size(a_perm,2));



figure; hold on;
cur_cv = 1; 
trueU = reshape(U(:,cur_cv),dur,size(U,1)/dur);
plot(nanmean(trueU,2),'r');
trueU = reshape(UU(:,cur_cv),dur,size(UU,1)/dur);
plot(nanmean(trueU,2),'b')
trueU = reshape(UUU(:,cur_cv),dur,size(UUU,1)/dur);
plot(nanmean(trueU,2),'g')


figure; hold on;
cur_cv = 1; 
trueU = reshape(V(:,cur_cv),dur,size(V,1)/dur);
plot(nanmean(trueU,2),'r');
trueU = reshape(VV(:,cur_cv),dur,size(VV,1)/dur);
plot(nanmean(trueU,2),'b')
trueU = reshape(VVV(:,cur_cv),dur,size(VVV,1)/dur);
plot(nanmean(trueU,2),'g')



a_temp = NaN(233,100);
for i = 1:100
  [a_full,b_full,~,U_full,V_full] = canoncorr(x_temp,y_temp); 
  a_temp(:,i) = a_full(:,1);
end


%% To do 
%choose two pairs of regions. 
%1) Find a motif in which MOs by self: 10
%2) motif that uses MOs and RSPd but not C: 5
%3) motif that uses All three: MOS, RSPd, VISp: 9 

%% In 10, test if any significant subspaces
%parse onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],10);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%mos and RSPD
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'RSPd'),:};
[~,~,~,~,r10_mos_rsp] = significantCVs(x,y,0.05,0);

%mos and visp
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'VISp'),:};
[~,~,~,~,r10_mos_vis] = significantCVs(x,y,0.05,0);

%vis and rsp
x = area_val{strcmp(area_label,'VISp'),:};
y = area_val{strcmp(area_label,'RSPd'),:};
[~,~,~,~,r10_vis_rsp] = significantCVs(x,y,0.05,0);

%% In 5
%parse onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],5);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%mos and RSPD
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'RSPd'),:};
[~,~,~,~,r5_mos_rsp] = significantCVs(x,y,0.05,0);

%mos and visp
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'VISp'),:};
[~,~,~,~,r5_mos_vis] = significantCVs(x,y,0.05,0);

%vis and rsp
x = area_val{strcmp(area_label,'VISp'),:};
y = area_val{strcmp(area_label,'RSPd'),:};
[~,~,~,~,r5_vis_rsp] = significantCVs(x,y,0.05,0);

%% In 9
%parse onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],9);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%mos and RSPD
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'RSPd'),:};
[~,~,~,~,r9_mos_rsp] = significantCVs(x,y,0.05,0);

%mos and visp
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'VISp'),:};
[~,~,~,~,r9_mos_vis] = significantCVs(x,y,0.05,0);

%vis and rsp
x = area_val{strcmp(area_label,'VISp'),:};
y = area_val{strcmp(area_label,'RSPd'),:};
[~,~,~,~,r9_vis_rsp] = significantCVs(x,y,0.05,0);


%% hypothesis; 
%MOS-RSP:
%strongest in 5 > 9 > 10
r5_mos_rsp
r9_mos_rsp
r10_mos_rsp

%MOS-VIS:
%strongest in 9 ... not really there in 5 or 10
r9_mos_vis
r5_mos_vis
r10_mos_vis

%VIS-RSP
%strongest in 9>5 ... not really there in 10
r9_vis_rsp
r5_vis_rsp
r10_vis_rsp

%% Next compare 3 and 7 MOs and BFD. Same spatial pattern but in 3 they occur in dequence and 7 synchrony
%parse onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],3);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%mos and RSPD
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'SSp-bfd'),:};
[a,b,U,V,r_3] = significantCVs(x,y,0.05,0);

%parse onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],7);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%mos and RSPD
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'SSp-bfd'),:};
[aa,bb,UU,VV,r_7] = significantCVs(x,y,0.05,0);

%Compare the map of weightings across neurons
figure; hold on; 
subplot(2,2,1); hold on; 
imagesc(a',[-0.3 0.3]); colormap(redgreencmap); colorbar
xlim([0.5,size(a,1)]);
title('MOs weights | Motif 3'); 

subplot(2,2,2); hold on;
imagesc(b',[-0.3 0.3]); colormap(redgreencmap); colorbar
xlim([0.5,size(b,1)]);
title('SS-BFD weights | Motif 7');

subplot(2,2,3); hold on;
imagesc(aa',[-0.3 0.3]); colormap(redgreencmap); colorbar
xlim([0.5,size(aa,1)]);
title('MOs weights | Motif 7');

subplot(2,2,4); hold on;
imagesc(bb',[-0.3 0.3]); colormap(redgreencmap); colorbar
xlim([0.5,size(bb,1)]);
title('SS-BFD weights | Motif 7');

sgtitle(sprintf('motif 3 CV: %0.2f,%0.2f,%0.2f motif 7: %0.2f ', r_3, r_7))

%write the leave-one-neuron out analysis

%% build a map for each neuron in both brain areas. 
%parse onsets
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],3);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%mos and RSPD
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'SSp-bfd'),:};
[delta_x, delta_y] = leaveOneOutCV(x,y,numel(r_3));


[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],7);

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');

%mos and RSPD
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'SSp-bfd'),:};
[delta_xx, delta_yy] = leaveOneOutCV(x,y,numel(r_7));

%Compare the map of weightings across neurons
close all; 
figure('position',[62 559 1179 420]); hold on; 
bar(delta_x*100,'EdgeAlpha',0,'FaceColor','r','FaceAlpha',0.5)
bar(delta_xx*100,'EdgeAlpha',0,'FaceColor','b','FaceAlpha',0.5)
title('MOs contributions | Motif 3 (r) vs 7 (b)'); 
ylabel('% Variance')
figure('position',[62 559 1179 420]); hold on; 
bar(delta_y*100,'EdgeAlpha',0,'FaceColor','r','FaceAlpha',0.5)
bar(delta_yy*100,'EdgeAlpha',0,'FaceColor','b','FaceAlpha',0.5)
title('SSp-bfd contributions | Motif 3 (r) vs 7 (b)'); 
ylabel('% Variance')

%% see the the relationship in activity of these subspaces during the opposite motif
%plot the MOs subspace (correct) during motif 3
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],3);
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'SSp-bfd'),:};
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
x = x-nanmean(x);
y = y-nanmean(y);

%true path 
figure; hold on; 
trueU = arrayfun(@(n) reshape(U(:,n),15,size(U,1)/15),1:size(U,2),'UniformOutput',0);
contraU = x*aa;
contraU = arrayfun(@(n) reshape(contraU(:,n),15,size(contraU,1)/15),1:size(contraU,2),'UniformOutput',0);
subplot(2,1,1); hold on;
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','r','LineWidth',1}), trueU)
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','b','LineWidth',1}), contraU)
title('MOs');
%BFD
trueV = arrayfun(@(n) reshape(V(:,n),15,size(V,1)/15),1:size(V,2),'UniformOutput',0);
contraV = y*bb;
contraV = arrayfun(@(n) reshape(contraV(:,n),15,size(contraV,1)/15),1:size(contraV,2),'UniformOutput',0);
subplot(2,1,2); hold on;
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','r','LineWidth',1}), trueV)
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','b','LineWidth',1}), contraV)
title('BFD');
sgtitle('Subpace activity during Motif 3')

%plot the MOs subspace (correct) during motif 7
[~,trig_st] = ParseByOnset([],st_norm,motif_onset,[1,15],7);
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'parent');
x = area_val{strcmp(area_label,'MOs'),:};
y = area_val{strcmp(area_label,'SSp-bfd'),:};
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
x = x-nanmean(x);
y = y-nanmean(y);

%true path 
figure; hold on; 
trueU = arrayfun(@(n) reshape(UU(:,n),15,size(UU,1)/15),1:size(UU,2),'UniformOutput',0);
contraU = x*a;
contraU = arrayfun(@(n) reshape(contraU(:,n),15,size(contraU,1)/15),1:size(contraU,2),'UniformOutput',0);
subplot(2,1,1); hold on;
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','r','LineWidth',1}), trueU)
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','b','LineWidth',1}), contraU)
title('MOs');
%BFD
trueV = arrayfun(@(n) reshape(VV(:,n),15,size(VV,1)/15),1:size(VV,2),'UniformOutput',0);
contraV = y*b;
contraV = arrayfun(@(n) reshape(contraV(:,n),15,size(contraV,1)/15),1:size(contraV,2),'UniformOutput',0);
subplot(2,1,2); hold on;
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','r','LineWidth',1}), trueV)
cellfun(@(x) shadedErrorBar(1:15,nanmean(x,2),sem(x,2),'lineprops',{'color','b','LineWidth',1}), contraV)
title('BFD');
sgtitle('Subpace activity during Motif 7')

%%

%%


























