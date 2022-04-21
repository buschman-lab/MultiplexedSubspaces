



%% Project the mean activity along either 
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace\Mouse 331 Recording 1regRRR_muaflag1_motif6.mat','st_norm','motif_onset','win','neu_area');

m = [5,14];
a = {'VIS','PRE','RSP'};

% m = [1,3];
% a = {'SSp-bfd','MOs','SS'};

[x,y] = loadFunc(neu_area,st_norm,motif_onset,win,m(1),[a(1,1),a(1,2)]); %subspace 1
[~,yy] = loadFunc(neu_area,st_norm,motif_onset,win,m(1),[a(1,1),a(1,3)]); %subspace 2

% [~,B,~] = RRR_simple(x,y);
% [~,BB,~] = RRR_simple(x,yy);

xx = nanmean(x,3)';
xt = preproX(x,x);
%% projection of motif average along either dimension
cur_t = 6;
a = xt*B(:,cur_t);
% a = normc(a);
b = xt*BB(:,cur_t);
% b = normc(b);

%summed variance 
figure; hold on; 
% plot(sum(nanvar(a,[],1)), sum(nanvar(b,[],1)),'Marker','x'); %activity along subspace 1 versus subspace 2
plot(a(:), b(:),'Marker','.','LineStyle','none'); %activity along subspace 1 versus subspace 2

a = qt*B(:,cur_t);
% a = normc(a);
b = qt*BB(:,cur_t);
% b = normc(b);

%summed variance 
% figure; hold on; 
% plot(sum(nanvar(a,[],1)), sum(nanvar(b,[],1)),'Marker','x'); %activity along subspace 1 versus subspace 2
plot(a(:), b(:),'Marker','.','LineStyle','none'); %activity along subspace 1 versus subspace 2


%% now do the other motif
a = {'VIS','PRE','RSP'};
% a = {'SSp-bfd','MOs','SS'};
[q,~] = loadFunc(neu_area,st_norm,motif_onset,win,m(2),[a(1,1),a(1,2)]); %subspace 1
[~,~] = loadFunc(neu_area,st_norm,motif_onset,win,m(2),[a(1,1),a(1,3)]); %subspace 2

% [~,B,~] = RRR_simple(x,y);
% [~,BB,~] = RRR_simple(x,yy);

qq = nanmean(q,3)';
qt = preproX(q,q);
%% 
a = qt*B(:,1:1);
% a = normc(a);
b = qt*BB(:,1:1);
% b = normc(b);

%summed variance 
% figure; hold on; 
% plot(sum(nanvar(a,[],1)), sum(nanvar(b,[],1)),'Marker','x'); %activity along subspace 1 versus subspace 2
plot(a(:), b(:),'Marker','.','LineStyle','none'); %activity along subspace 1 versus subspace 2


%%

figure; hold on 
for i = 1:10
a = xt*B(:,i);
% a = normc(a);
b = xt*BB(:,i);
% b = normc(b);
plot(sum(nanvar(a,[],1)), sum(nanvar(b,[],1)),'Marker','x','color','k'); %activity along subspace 1 versus subspace 2
end
% 
set(gca,'XScale','log','YScale','log')
xvals = get(gca,'xlim');
yvals = get(gca,'ylim');
plot(xvals,yvals)




%% compare motifs
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace\Mouse 331 Recording 1regRRR_muaflag1_motif6.mat','st_norm','motif_onset','win','neu_area');

m = [5,14];
a = {'VIS','PRE','RSP'};

% m = [1,3];
% a = {'SSp-bfd','MOs','SS'};


%% load the activity
[x,y] = loadFunc(neu_area,st_norm,motif_onset,win,m(1),a); %motif 1
[xx,yy] = loadFunc(neu_area,st_norm,motif_onset,win,m(2),a); %motif 2
% Get same number of trials in both motifs
% n=min(size(x,3),size(xx,3));
% idx = randperm(size(x,3),n); %so same random trials in x and y
% x = x(:,:,idx);
% y = y(:,:,idx);
% idx = randperm(size(xx,3),n); %so same random trials in x and y
% xx = xx(:,:,idx);
% yy = yy(:,:,idx);

%normalize the variances 
%% hold out 30% of trials and combined
rng('default');
n = size(x,3);
testidx = randperm(n,floor(0.2*n));
z = numel(testidx);
xt = nanmean(x,3)';
% yt = nanmean(abs(x-nanmean(x,3)),3)';
yt = preproX(x(:,:,testidx),x(:,:,testidx));
trainidx = ismember(1:n, testidx)==0;
[~,B,V] = RRR_simple(x(:,:,trainidx),y(:,:,trainidx));
% 
n = size(xx,3);
testidx = randperm(n,floor(0.2*n));
xxt = nanmean(xx,3)';
yyt = preproX(xx(:,:,testidx),xx(:,:,testidx));
% yyt = nanmean(abs(xx-nanmean(xx,3)),3)';
trainidx = ismember(1:n, testidx)==0;
[~,BB,VV] = RRR_simple(xx(:,:,trainidx),yy(:,:,trainidx));
zz = numel(testidx);


%% test if the variance along the subspace dims 









%% go through withheld trials in both and project along 
a = reshape(yt*B,z,12,size(B,2));
b = reshape(yt*BB,z,12,size(BB,2));
n = min(zz,z);
a = squeeze(nanmean(a(1:n,:,:),2));
b = squeeze(nanmean(b(1:n,:,:),2));

plot(a(:),b(:),'.')


%%
a = yt*B(:,1)


plot(a(:),b(:),'.')
corr(a(:),b(:))




%% Is the average pattern of activity in motif n and m the same in the source area
a = nanmax(xt,[],1);
b = nanmax(xxt,[],1);
plot([0,1.8],[0,1.8],'color','k')
figure; hold on;
plot(a(:),b(:),'.')
corr(a(:),b(:)) 
xlabel('motif n trial average ');
ylabel('motif m trial average ');


%% Project 
a = nanmean(xt,1)*B(:,1:10);
b = nanmean(xt,1)*BB(:,1:10);
c = nanmean(xxt,1)*B(:,1:10);
d = nanmean(xxt,1)*BB(:,1:10);
figure; hold on;
plot(a(:),b(:),'.','color','k')
plot(c(:),d(:),'.','color','r')
corr(a(:),b(:)) 



%% we can also do this per trial
a = nanmean(abs(yt),1);
b = nanmean(abs(yyt),1);
figure; hold on;
plot(a(:),b(:),'.')
corr(a(:),b(:))
xlabel('motif n noise');
ylabel('motif m noise');

%% Is our motifs activity alligned with our average noise along our subspace
a = xt*B;
b = yt*B;

figure; hold on;
for i = 1:10
    plot(a(:,i),b(:,i),'.'); 
end
xlabel('motif n subspace n ');
ylabel('motif n subspace n noise ');

%% If we take the other motifs activity and project it along this subspace
c = xxt*B;
figure; hold on;
for i = 1:10
    plot(a(:,i),c(:,i),'.'); 
end
xlabel('motif n subspace n ');
ylabel('motif m subspace n ');


%% If we take the other motifs noise and project it along this subspace
d = yyt*B;
figure; hold on;
for i = 1:10
    plot(b(:,i),d(:,i),'.'); 
end
xlabel('motif n subspace n noise');
ylabel('motif m subspace n noise');

%% What if we take Motif N subspace N and motif N subspace M
a = xt*B;
b = xt*BB;
figure; hold on;
for i = 1:10
    plot(a(:,i),b(:,i),'.'); 
end
xlabel('motif n subspace n');
ylabel('motif n subspace m');

%% What if we take Motif N subspace N and motif M subspace M
a = xt*B;
b = xxt*BB;
figure; hold on;
for i = 1:10
    plot(a(:,i),b(:,i),'.'); 
end
xlabel('motif n subspace n');
ylabel('motif m subspace m');


%% What if we take Motif N subspace N and motif N subspace M
a = xt*B(:,1:10);
b = xt*BB(:,1:10);

figure; hold on;
for i = 1:10
    plot(a(:,i),b(:,i),'.'); 
end
xlabel('motif n subspace n');
ylabel('motif n subspace m');

%%
a = nanmean(xt)*B(:,1:10);
b = nanmean(xt)*BB(:,1:10);

c = nanmean(xxt)*B(:,1:10);
d = nanmean(xxt)*BB(:,1:10);

figure; hold on;
plot(a,b,'o'); 
plot(c,d,'o'); 
xlabel('motif n subspace n');
ylabel('motif n subspace m');









%% get average activity per trial
rng('default');
n = size(x,3);
testidx = randperm(n,floor(0.2*n));
a = squeeze(nanmean(x(:,:,testidx),[2,3]));
n = size(xx,3);
testidx = randperm(n,floor(0.2*n));
b = squeeze(nanmean(xx(:,:,testidx),[2,3]));
figure; hold on;
plot(a(:),b(:),'.')
corr(a(:),b(:))




%% now reorder
[rho,idx] = arrayfun(@(n) max(abs(corr(B(:,n),BB))), 1:10);
B_re = B(:,idx);
B = B(:,1:10);


%%
a = xtest*B(:,4);
b = xtest*B_re(:,4);
c = xxt*B(:,1);
d = xxt*B_re(:,1);

figure; hold on; 
plot(a,b,'.')




% [~,B,V] = RRR_simple(xtrain(1:137,:,:),ytrain(1:size(ytrain,1)/2,:,:));
% [~,BB,VV] = RRR_simple(xtrain(138:end,:,:),ytrain((1+(size(ytrain,1)/2)):end,:,:));


%%
[x,y] = preproX(xtest(1:137,:,:),ytest(1:137,:,:));
%%
% y = x*B(:,1)*V(:,1)';
% yhat = x*BB(:,1)*VV(:,1)';

for i = 1:10
    y = x*B(:,i)*V(:,i)';
    yhat = x*BB(:,i)*VV(:,i)';
    pev(i) = 1-NormalizedSquaredError(y,yhat);
end

%% project the withheld
%split our betas by target 
imagesc(rrr_B(:,1:10));
ndim = 2;
n = size(rrr_B,1)/2;
B = rrr_B(1:n,1:ndim);
BB = rrr_B((n+1):end,1:ndim);
cur_t = 1; 
[x,y] = preproX(xtest,ytest);
%%
xdim1 = x*B(:,1:ndim);
xdim2 = x*B(:,1:ndim);






%%

function [x,y] = loadFunc(neu_area,st_norm,motif_onset,win,motif,area_name)
    [~,trig_st] = ParseByOnset([],st_norm,motif_onset,win,motif);

    %parse activity per parent region 
    [area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'general');

    %remove the rare edge case where a motif begins at the start (no baseline)
    area_val = RemoveEdgeTrials(area_val);

    %clean up areas %third input is the min # of spikes to keep area
    [area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 
    
    x = area_val{strcmp(area_label,area_name{1})};
%     y = cat(1,area_val{strcmp(area_label,area_name{2})},area_val{strcmp(area_label,area_name{3})});
    y = area_val{strcmp(area_label,area_name{2})};

    %normalize to baseline
    x = normalizeToBaseline(x,[1:2],'mean');
    y = normalizeToBaseline(y,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    y = y(:,3:end,:);
end

function [x,y] = preproX(x,y)
%subtract the psth
x = x-nanmean(x,3);
y = y-nanmean(y,3);

%concatentate across trials and pca
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

end
