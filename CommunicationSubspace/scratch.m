data = LoadSubspaceData('paired');







%% Plot the PCA projection 
cur_rec = 6;
cur_motif = [5,6];
area  = {'VIS','PRE','RSP'};
% cur_motif = [1,3];
% area = {'SSp-bfd','MOs','SS'};
RunPCA(data,cur_rec,cur_motif,area(1));


%% How does this relate to the subspaces? 
%choose the dimensions with greatest difference between the two 
%the do pca of one population and 




%%
% cur_motif = [1,3];
% area = {'SSp-bfd','MOs','SS'};
cur_motif = [5,14];
% area = {'SSp-bfd','MOs','SS'};
area  = {'VIS','PRE','RSP'};
cur_rec = [1,2,4,5,6];
a = NaN(4,2);
b = NaN(4,2);
for j = 1:2
    for i = 1:numel(cur_rec)
        x = loadFunc(data,cur_rec(i),cur_motif(j),area{1},0); 
        [B,~] = loadCoef(data,cur_rec(i),cur_motif(j),{area{1},area{2}});
        [BB,~] = loadCoef(data,cur_rec(i),cur_motif(j),{area{1},area{3}});
        %get the subspace variance projection of the mean activity
        xx = nanmean(x,3)';
        a(i,j) = sum(nanvar(xx*B));
        b(i,j) = sum(nanvar(xx*BB));
    end
end

%%
figure; hold on; 
plot(a(:,1),b(:,1),'x','color','k')
plot(a(:,2),b(:,2),'x','color','g')

xvals = get(gca,'xlim');
yvals = get(gca,'ylim');
plot([0,max([xvals,yvals])],[0,max([xvals,yvals])])

%% grab an example recording and projection
area  = {'VIS','PRE','RSP'};
x = loadFunc(data,1,14,area{1},0); 
[B,~] = loadCoef(data,1,14,{area{1},area{2}});
[BB,~] = loadCoef(data,1,14,{area{1},area{3}});
%get the subspace variance projection of the mean activity
xx = nanmean(x,3)';
a = xx*B;
b = xx*BB;

figure; hold on; 
subplot(1,2,1);
plot(a); title('projecting to PRE');
xlabel('trial timepoints');
subplot(1,2,2);
plot(b); title('projecting to RSP');
xlabel('trial timepoints');
sgtitle('motif 14');

%% 
area  = {'VIS','PRE','RSP'};
x = loadFunc(data,1,5,area{1},0); 
[B,~] = loadCoef(data,1,5,{area{1},area{2}});
[BB,~] = loadCoef(data,1,5,{area{1},area{3}});
%get the subspace variance projection of the mean activity
xx = nanmean(x,3)';
a = xx*B;
b = xx*BB;

figure; hold on; 
subplot(1,2,1);
plot(a); title('projecting to PRE');
xlabel('trial timepoints');
subplot(1,2,2);
plot(b); title('projecting to RSP');
sgtitle('motif 5');
xlabel('trial timepoints');
%% 














%%
function [B,V,lambda] = loadCoef(data,cur_rec,cur_motif,area_name) 
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    paired_areas = data{cur_rec}(cur_motif).paired_areas;
    
    a = find(strcmp(area_label,area_name{1}));
    b = find(strcmp(area_label,area_name{2}));
    idx = find(ismember(paired_areas,[a,b],'rows')); 
    
    B = data{cur_rec}(cur_motif).rrr_B{idx};
    V = data{cur_rec}(cur_motif).rrr_V{idx};
    
    %get the originally used lambda
    x = area_val{strcmp(area_label,area_name{1})};
    if size(x,3)>1
        x = x-nanmean(x,3);
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';        
    end
    dMaxShrink = .5:.01:1;
    lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);
    cvl_ridge = data{cur_rec}(cur_motif).cvl_ridge{idx};
    [~,idx] = bestLambda(cvl_ridge);
    lambda = lambda(idx);
    
end


function x = loadFunc(data,cur_rec,cur_motif,area_name,flag)
    if nargin <5; flag = 1; end
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    x = area_val{strcmp(area_label,area_name)};

    %normalize to baseline
%     x = normalizeToBaseline(x,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    
    if flag ==1 
        %subtract the psth
        x = x-nanmean(x,3);

        %concatentate across trials and pca
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    elseif flag == 2
        %concatentate across trials and pca
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    end
end


function RunPCA(data,cur_rec,cur_motif,area)
%get all motifs
rng('default');
z = arrayfun(@(n) loadFunc(data,cur_rec,n,area{1},0),1:14,'UniformOutput',0); 
n = min(cellfun(@(x) size(x,3),z));
testidx = randperm(n,floor(0.3*n));
trainidx = ismember(1:n, testidx)==0;
a = z{cur_motif(1)}(:,:,testidx);
% b = cellfun(@(x) x(:,:,testidx), z(ismember(1:14,cur_motif(1))==0),'UniformOutput',0);
% b = cat(3,b{:});
b = z{cur_motif(2)}(:,:,testidx);
z = cellfun(@(x) x(:,:,trainidx),z,'UniformOutput',0);
% a = loadFunc(data,cur_rec,cur_motif(1),area{1},0); 
% b = loadFunc(data,cur_rec,cur_motif(2),area{1},0); 

%get subset of data from each motif
% % n = min(size(a,3),size(b,3));
% idx = randperm(size(a,3),n); %so same random trials in x and y
% a = a(:,:,idx);
% idx = randperm(size(b,3),n); %so same random trials in x and y
% b = b(:,:,idx);
% testidx = randperm(n,floor(0.3*n));
% trainidx = ismember(1:n, testidx)==0;
%% training data
% x = cat(3,a(:,:,trainidx),b(:,:,trainidx));
% x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
% x = x-nanmean(x); %mean center
x = cat(3,z{:});
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
x = x-nanmean(x); %mean center
[coef,score] = pca(x);

%project withheld trials
y = cell(numel(size(a,3)),1);
for i = 1:size(a,3)
    x = a(:,:,i)';
    x = x-nanmean(x);
    y{i} = x*coef(:,1:4);
end
yy = cell(numel(size(b,3)),1);
for i = 1:size(b,3)
    x = b(:,:,i)';
    x = x-nanmean(x);
    yy{i} = x*coef(:,1:4);    
end

%% figure plot

a = permute(cat(3,y{:}),[2,1,3]);
aa = nanmean(a,3)';
a = reshape(a,size(a,1),size(a,2)*size(a,3))';
b = permute(cat(3,yy{:}),[2,1,3]);
bb = nanmean(b,3)';
b = reshape(b,size(b,1),size(b,2)*size(b,3))';

close all
figure('units','normalized','position',[0.1 0.1 0.8 0.4]); hold on;  
subplot(1,4,1); hold on 
xval = cat(1,a(:,1),b(:,1));
plot(a(:,1),a(:,2),'marker','.','linestyle','none','color','r','markersize',1); AddLSline(a(:,1),a(:,2),xval,[0.8 0.1 0.1],1.5);
plot(b(:,1),b(:,2),'marker','.','linestyle','none','color','b','markersize',1); AddLSline(b(:,1),b(:,2),xval,[0.1 0.1 0.8],1.5);
xlabel('dim 1'); ylabel('dim 2');
subplot(1,4,2); hold on
xval = cat(1,a(:,1),b(:,1));
plot(a(:,1),a(:,3),'marker','.','linestyle','none','color','r','markersize',1); AddLSline(a(:,1),a(:,3),xval,[0.8 0.1 0.1],1.5); 
plot(b(:,1),b(:,3),'marker','.','linestyle','none','color','b','markersize',1); AddLSline(b(:,1),b(:,3),xval,[0.1 0.1 0.8],1.5);
xlabel('dim 1'); ylabel('dim 3');
subplot(1,4,3); hold on
xval = cat(1,a(:,2),b(:,2));
AddLSline(a(:,2),a(:,3),xval,[0.8 0.1 0.1],1.5); 
AddLSline(b(:,2),b(:,3),xval,[0.1 0.1 0.8],1.5);
plot(a(:,2),a(:,3),'marker','.','linestyle','none','color','r','markersize',1); 
plot(b(:,2),b(:,3),'marker','.','linestyle','none','color','b','markersize',1); 
xlabel('dim 2'); ylabel('dim 3');
subplot(1,4,4); hold on
xval = cat(1,a(:,3),b(:,3));
AddLSline(a(:,3),a(:,4),xval,[0.8 0.1 0.1],1.5); 
AddLSline(b(:,3),b(:,4),xval,[0.1 0.1 0.8],1.5);
plot(a(:,3),a(:,4),'marker','.','linestyle','none','color','r','markersize',1); 
plot(b(:,3),b(:,4),'marker','.','linestyle','none','color','b','markersize',1); 
xlabel('dim 3'); ylabel('dim 4');

yy = cat(3,yy{:});
y = cat(3,y{:});
% figure('units','normalized','position',[0.1 0.1 0.8 0.4]); hold on;  
% subplot(1,3,1); hold on
% cellfun(@(x) plot3(1:12,x(:,1),x(:,2),'marker','none','linestyle','-','color','r','markersize',5), {nanmean(y,3)});
% cellfun(@(x) plot3(1:12,x(:,1),x(:,2),'marker','none','linestyle','-','color','b','markersize',5), {nanmean(yy,3)});
% set(gca,'CameraPosition',[-40.5762 -3.3029   0.4629]);
% xlabel('time'); ylabel('dim 1'); zlabel('dim 2');
% axis square
% subplot(1,3,2); hold on
% cellfun(@(x) plot3(1:12,x(:,1),x(:,3),'marker','none','linestyle','-','color','r','markersize',5), {nanmean(y,3)});
% cellfun(@(x) plot3(1:12,x(:,1),x(:,3),'marker','none','linestyle','-','color','b','markersize',5), {nanmean(yy,3)});
% set(gca,'CameraPosition',[-40.5762 -3.3029   0.4629]);
% xlabel('time'); ylabel('dim 1'); zlabel('dim 3');
% axis square
% subplot(1,3,3); hold on
% cellfun(@(x) plot3(1:12,x(:,2),x(:,3),'marker','none','linestyle','-','color','r','markersize',5), {nanmean(y,3)});
% cellfun(@(x) plot3(1:12,x(:,2),x(:,3),'marker','none','linestyle','-','color','b','markersize',5), {nanmean(yy,3)});
% set(gca,'CameraPosition',[-40.5762 -3.3029   0.4629]);
% xlabel('time'); ylabel('dim 2'); zlabel('dim 3');
% axis square

figure; hold on; 
subplot(1,3,1); hold on
shadedErrorBar(1:12,nanmean(y(:,1,:),3),sem(y(:,1,:),3),'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
shadedErrorBar(1:12,nanmean(yy(:,1,:),3),sem(yy(:,1,:),3),'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});
subplot(1,3,2); hold on
shadedErrorBar(1:12,nanmean(y(:,2,:),3),sem(y(:,2,:),3),'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
shadedErrorBar(1:12,nanmean(yy(:,2,:),3),sem(yy(:,2,:),3),'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});
subplot(1,3,3); hold on
shadedErrorBar(1:12,nanmean(y(:,3,:),3),sem(y(:,3,:),3),'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
shadedErrorBar(1:12,nanmean(yy(:,3,:),3),sem(yy(:,3,:),3),'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});

figure; hold on; 
subplot(1,3,1); hold on
% y = cat(3,y{:});
shadedErrorBar(nanmean(y(:,1,:),3),nanmean(y(:,2,:),3),sem(y(:,2,:),3),'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
% yy = cat(3,yy{:});
shadedErrorBar(nanmean(yy(:,1,:),3),nanmean(yy(:,2,:),3),sem(yy(:,2,:),3),'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});
xlabel('dim 1'); ylabel('dim 2');
axis square
subplot(1,3,2); hold on
shadedErrorBar(nanmean(y(:,2,:),3),nanmean(y(:,3,:),3),sem(y(:,3,:),3),'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
shadedErrorBar(nanmean(yy(:,2,:),3),nanmean(yy(:,3,:),3),sem(yy(:,3,:),3),'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});
xlabel('dim 2'); ylabel('dim 3');
axis square
subplot(1,3,3); hold on
shadedErrorBar(nanmean(y(:,1,:),3),nanmean(y(:,3,:),3),sem(y(:,3,:),3),'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
shadedErrorBar(nanmean(yy(:,1,:),3),nanmean(yy(:,3,:),3),sem(yy(:,3,:),3),'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});
xlabel('dim 1'); ylabel('dim 3');
axis square



end
















