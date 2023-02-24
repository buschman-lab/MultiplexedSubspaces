function x = Plot_MotifAlignment(data)
%Camden MacDowell
%data is paired data | data = LoadSubspaceData('paired');
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Rotations';

%%
close all
cur_rec = 6;
cur_motif = [5,6];
area  = {'VIS','PRE','SSp-bfd'}; %important | this should always go Starting Point --> Higher in Motif A --> Higher in Motif B; 
[x,xx] = MainFunc(data,cur_rec,cur_motif,area);
% saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('motif5_6'),savedir,0); close all
% saveCurFigs(gcf,{'-dsvg','-dpng'},sprintf('motif5_6_remake'),savedir,0); close all

%%

close all;
rng('default')
cur_rec = 6;
cur_motif = [8,6];
area  = {'SSp-bfd','SS','RSP'};
[y,yy] = MainFunc(data,cur_rec,cur_motif,area);
% saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('motif8_6_subspaces'),savedir,0); close all
% saveCurFigs(gcf,{'-dsvg','-dpng'},sprintf('motif8_6_remake'),savedir,0); close all
%%
cur_motif = [3,14];
area  = {'MOs','SSp-bfd','THAL'};
[z,zz] = MainFunc(data,cur_rec,cur_motif,area);
% saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('motif1_13'),savedir,0); close all
% saveCurFigs(gcf,{'-dsvg','-dpng'},sprintf('motif3_14_remake'),savedir,0); close all
%% 
close all
cur_motif = [3,9];cc
area  = {'MOs','SS','VIS'}; 
[q,qq] = MainFunc(data,cur_rec,cur_motif,area);
% saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('motif3_9'),savedir,0); close all
% saveCurFigs(gcf,{'-dsvg','-dpng'},sprintf('motif3_9_remake'),savedir,0); close all
%% 
cur_motif = [14,6];
area  = {'VIS','HIPP','SSp-bfd'};
[r,rr] = MainFunc(data,cur_rec,cur_motif,area);
% saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('motif5_9'),savedir,0); close all
% saveCurFigs(gcf,{'-dsvg','-dpng'},sprintf('motif14_6_remake'),savedir,0); close all
%%
cur_motif = [4,11];
area  = {'MOs','PRE','SSp-bfd'}; %could also do HIPP
[w,ww] = MainFunc(data,cur_rec,cur_motif,area);
% saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('motif4_11'),savedir,0); close all
% saveCurFigs(gcf,{'-dsvg','-dpng'},sprintf('motif4_11_remake'),savedir,0); close all

%%  Plot the combined statistics
theta_all = {x(:),y(:),z(:),q(:),r(:),w(:)};
theta_all = [{cat(1,theta_all{:})},theta_all];
xpos=[1,2.5,3.5,4.5,5.5,6.5,7.5];

theta_all = theta_all([1,2,6,3,7,5,4]);
plotCombinedSubspaceAngles(theta_all,xpos)
%%
theta_all = {x(:),y(:),z(:),q(:),r(:),w(:)};
plotSubspaceAngles(cat(1,theta_all{:}))
%%
saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('CombinedStatistics_new'),savedir,0); close all

%% Run across all pairs of motifs and all pairs of recordings 
%make sure recording actually fits the criteria that 1) they significantly
%differ in at least one region (ttest) and that the direction is correct
%(otherwise the angles will be reverse) such that A>B for the first 'target' and
%B>2 for the second 'target'

%all recordings
%all triads of areas
a = nchoosek(1:8,3);
%get all orderings
a_all = NaN((size(a,1)*3),3);
COUNT = 1;
for i = 1:size(areas,1)    
    for j = 1:3
        a_all(COUNT,:) = a(i,circshift([1,2,3],j-1));
        COUNT = COUNT+1;   
    end
end
%all pairs of motifs
m_all = nchoosek(1:14,2);
m_all = cat(1,m_all,fliplr(m_all));

area_label = {'HIPP','MOs','PRE','RSP','SS','SSp-bfd','THAL','VIS'};
combined_data = [];
COUNT = 1; 
for cur_rec = 6    
    for a = 1:size(a_all,1)
        a
        for m = 1:size(m_all,1)
           %check criteria
           [crit,sig] = CompareActivity(data,cur_rec,m_all(m,:),area_label(a_all(a,:)),0);
           %is A>B and then B>A and at least one sig diff?
           if crit(1)==1 && crit(2) == 0 && sum(sig)>0 %then check the angles across recordings
              x = GetSubspaceAngleForRecordings(data,m_all(m,:),area_label(a_all(a,:)),cur_rec);
              x = cat(1,x{:});
              temp = x(:,1)-x(:,2);
              combined_data(COUNT).d = temp;
              combined_data(COUNT).rec = cur_rec;
              combined_data(COUNT).area = a_all(a,:);
              combined_data(COUNT).pair = m_all(m,:);
              COUNT = COUNT+1;             
           end           
        end
    end
end



end

%% subroutines
function plotCombinedSubspaceAngles(xx,xpos)
fp = fig_params_cortdynamics;
figure; hold on;
% combined stats
for i = 1:numel(xx)    
    x = cat(1,xx{i}{:});
    if size(x,2)>1
        theta = x(:,1)-x(:,2);
    else
        theta=x; 
    end
    %get a bootstrapped distribution
    [ci,thetaboot] = bootci(1000,@circ_mean, deg2rad(theta));
    ci = rad2deg(ci);
    mu = rad2deg(circ_mean(thetaboot));
    err = rad2deg(circ_std(thetaboot));
    thetaboot = rad2deg(thetaboot);
    p = sum(thetaboot>=0)/numel(thetaboot);
    text(xpos(i),7,sprintf('%0.2f +/- %0.2f \n[%0.2f %0.2f] \n p=%0.4f',mu,err,ci(1),ci(2),p),'FontSize',5);
    if i==1
        CompareViolins(thetaboot',fp,'plotspread',0,'divfactor',0.2,'xpos',xpos(i),'col',{fp.c_ff},'distWidth',1.4);    
    else
        CompareViolins(thetaboot',fp,'plotspread',0,'divfactor',0.5,'xpos',xpos(i),'col',{[0.25 0.25 0.25]});
    end
end
camroll(-90)
set(gca,'xtick',xpos)
xlim([0 xpos(end)+0.75])
plot([0 xpos(end)+1],[0 0],'linestyle','-','color',[0.3 0.3 0.3],'linewidth',1.5)
yvals = get(gca,'ylim');
% ylim([floor(yvals(1)),max(ceil(yvals(2)),5)]);
ylim([-10 5]);
ylabel('\Delta theta');
set(gca, 'YGrid', 'off', 'XGrid', 'on','GridAlpha',0.15);
fp.FormatAxes(gca); box on; 
fp.FigureSizing(gcf,[3 3 4 5.75],[10 10 14 10])
end



function [x,xx] = MainFunc(data,cur_rec,cur_motif,area)
a = cell(1,3);
b = cell(1,3);
a_all = cell(1,3);
b_all = cell(1,3);
theta = cell(1,3);
thetacon = cell(1,3);
for i = 1:3
    [~,~,y,yy,aa,bb,yyy] = RunPCA(data,cur_rec,cur_motif,area(i)); %y is the Motif A activity along Motif A PC1, YY is Motif B activity along motif A PC1, yyy is Motif B activity along Motif B PC1 (which is what we use for comparison later)
    PlotExampleSubspaceRelationships(data,cur_rec,cur_motif,[area(i),area(ismember(1:3,i)==0)],y,yy,aa,bb,yyy)
%     PlotProjection(y,yy,cur_motif,area(i),cur_rec,[1,2]);
    a{i} = nanmean(y(:,1:3,:),3);
    b{i} = nanmean(yy(:,1:3,:),3);
    a_all{i} = squeeze(y(:,1,:));
    b_all{i} = squeeze(yy(:,1,:));
    %compare to subspace betas
    [theta{i},thetacon{i}] = CompareToSubspace(data,cur_rec,cur_motif,[area(i),area(ismember(1:3,i)==0)],1000,i); %
end
% compare the joint projections 
PlotJointProjections(a,b,cur_motif,area,cur_rec);
% get bootstrapped angle between projection
JointProjectionAngle(a_all,b_all); %
% compare these two motifs trialwise activity in each region
CompareActivity(data,cur_rec,cur_motif,area,1);  
% plot the anges
CompareAngles(theta,thetacon,area)
% angle between joint distribution
% thetaActivity = NaN(2,6);
% for cur_rec = 1:6
%     try
%     [thetaActivity(1,cur_rec), thetaActivity(2,cur_rec)] = GetActivityAngleAcross(data,cur_motif,area,cur_rec);    
%     catch
%     end
% end
% 
% %get the subspace 
% for cur_rec = 1:6
%     try 
%     catch 
%     end
% end
% plot angles subspace--> activity across recordings pairwise
x = GetSubspaceAngleAcrossRecordings(data,cur_motif,area);

plotSubspaceAngles(x)
%average across recordings
xx = cell(size(x,1),1);
for j = 1:size(x,1)
   temp = cat(1,x{j,:});
   temp = temp(:,1)-temp(:,2);
   temp = rad2deg(circ_mean(deg2rad(temp)));
   xx{j} = temp;
end


end %function end




function theta = JointProjectionAngle(a,b)
fp = fig_params_cortdynamics;
rng('default');
theta = NaN(1,1000);
for i = 1:1000 %bootstraps
   x = cat(3,a{:});
   y = cat(3,b{:});
   idx = datasample(1:size(x,2),size(x,2)); 
   x = squeeze(nanmean(x(:,idx,:),2));   
   idx = datasample(1:size(y,2),size(y,2));
   y = squeeze(nanmean(y(:,idx,:),2));
   theta(i) = rad2deg(subspace(x,y));
end

figure;
[mu,err] = AnglePlot(theta,[0.2 0.2 0.2]);
fp.FigureSizing(gcf,[3 3 2 2],[10 10 14 10])
title(sprintf('%0.2f +/- %0.2f',mu,err),'FontWeight','normal')
end 

function [rho] = PlotExampleSubspaceRelationships(data,cur_rec,cur_motif,area,y,yy,a,b,yyy,verbose)
if nargin <10; verbose = 1; end
fp = fig_params_cortdynamics;
%relative to the subspaces of motif A 
B = loadCoef(data,cur_rec,cur_motif(1),area([1,2])); 
BB = loadCoef(data,cur_rec,cur_motif(2),area([1,2])); 

%match the direction of the Betas
if corr(B(:,1),BB(:,1))<corr(B(:,1),-1*BB(:,1))
    BB = -1*BB;
end

rng('default');
%get the angle of rotation of the subspace weights
x = reshape(a,[size(a,1),size(a,2)*size(a,3)])';
x = x-nanmean(x); %mean center
coef = pca(x);
subtheta =[];
subtheta(1) = rad2deg(subspace(coef(:,1),B(:,1)));
x = reshape(b,[size(b,1),size(b,2)*size(b,3)])';
x = x-nanmean(x); %mean center
coef = pca(x);
subtheta(2) = rad2deg(subspace(coef(:,1),B(:,1)));
subtheta(3) = rad2deg(subspace(coef(:,1),BB(:,1)));

%get the average projection into subspace
proj1 = [];
x = nanmean(a,3)';
proj1(:,1) = nanmean(squeeze(y(:,1,:)),2); %average projection in the PC1 (RED)
proj1(:,2) = x*B(:,1);

%do for motif B on it's own PC projected along the subspace from Motif A
%(BLUE) - this is what we use later
proj2 = [];
x = nanmean(b,3)';
proj2(:,1) = nanmean(squeeze(yyy(:,1,:)),2);
proj2(:,2) = x*B(:,1);

%do for motif B on it's own PC projected along the subspace from Motif B (GREEN)
proj3 = [];
x = nanmean(b,3)';
proj3(:,1) = nanmean(squeeze(yyy(:,1,:)),2); 
proj3(:,2) = x*BB(:,1);

% % %standardize the subspace to zero
% proj1(:,2) = proj1(:,2)-min(proj1(:,2));
% proj2(:,2) = proj2(:,2)-min(proj2(:,2));
% proj3(:,2) = proj3(:,2)-min(proj3(:,2));
% % 
% %standardize the subspace to zero
% proj1(:,1) = proj1(:,1)-min(proj1(:,1));
% proj2(:,1) = proj2(:,1)-min(proj2(:,1));
% proj3(:,1) = proj3(:,1)-min(proj3(:,1));


if verbose        
    figure; hold on;
    %plot the line of best fit
    m = [];
    [m(1),xint] = AddLSline(proj1(:,1),proj1(:,2),proj1(:,1),[0.8 0.1 0.1],1.5);
    m(2) = AddLSline(proj2(:,1),proj2(:,2),proj2(:,1),[0.1 0.1 0.8],1.5);
    m(3) = AddLSline(proj3(:,1),proj3(:,2),proj3(:,1),[0.1 0.8 0.1],1.5);
    plot(proj1(:,1),proj1(:,2),'o','color',[0.8 0.1 0.1],'markersize',fp.markersizesmall,'markerfacecolor',[0.8 0.1 0.1],'linestyle','none','linewidth',1.5);
    plot(proj2(:,1),proj2(:,2),'o','color',[0.1 0.1 0.8],'markersize',fp.markersizesmall,'markerfacecolor',[0.1 0.1 0.8],'linestyle','none','linewidth',1.5);
    plot(proj3(:,1),proj3(:,2),'o','color',[0.1 0.8 0.1],'markersize',fp.markersizesmall,'markerfacecolor',[0.1 0.8 0.1],'linestyle','none','linewidth',1.5);
    %add a line at 45degrees 
    fx = @(x) 1*x+xint;
    unity = fx(proj1(:,1));
    plot(proj1(:,1),unity,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle','-');
    fp.FormatAxes(gca)
    fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
    
%     xvals = get(gca,'xlim');
%     yvals = get(gca,'ylim');
%     set(gca,'xlim',[-0.1*xvals(2),xvals(2)]);
%     set(gca,'ylim',[-0.1*yvals(2),yvals(2)]);
    
    %add lines showing the distribution along each axis
    xvals = get(gca,'xlim');
    plot([xvals(1),xvals(1)],[min(proj1(:,2)),max(proj1(:,2))],'linewidth',3,'color',[0.8 0.1 0.1])
    plot([xvals(1),xvals(1)],[min(proj2(:,2)),max(proj2(:,2))],'linewidth',3,'color',[0.1 0.1 0.8])
    plot([xvals(1),xvals(1)],[min(proj3(:,2)),max(proj3(:,2))],'linewidth',3,'color',[0.1 0.8 0.1])
    yvals = get(gca,'ylim');
    plot([min(proj1(:,1)),max(proj1(:,1))],[yvals(1),yvals(1)],'linewidth',3,'color',[0.8 0.1 0.1])
    plot([min(proj2(:,1)),max(proj2(:,1))],[yvals(1),yvals(1)],'linewidth',3,'color',[0.1 0.1 0.8])
    plot([min(proj3(:,1)),max(proj3(:,1))],[yvals(1),yvals(1)],'linewidth',3,'color',[0.1 0.8 0.1])
    
    %get the angle between the lines
    t = nchoosek(1:numel(m),2);
    theta = NaN(1,size(t,1));
    for i = 1:numel(m)
        aa = m(t(i,1));
        bb = m(t(i,2));
        theta(i) = abs(rad2deg(atan(aa))-rad2deg(atan(bb))); %abs(rad2deg(atan( (aa-bb)/(1+aa*bb))));
    end
    xlabel(sprintf('Population Activity (PC 1)\n [%s]',area{1}));
    ylabel(sprintf('Subspace Activity (Dim 1)\n [%s-%s]',area{1},area{2}));
    title(sprintf('\ntheta = r-->b %0.2f \nr-->g %0.2f b-->g %0.2f',theta),'fontweight','normal');                
    
    figure;
    AnglePlot(rad2deg(atan(m(1))),[0.8 0.1 0.1]) 
    AnglePlot(rad2deg(atan(m(2))),[0.1 0.1 0.8]) 
    AnglePlot(rad2deg(atan(m(3))),[0.1 0.8 0.1]) 
    AnglePlot(45,[0.1 0.1 0.1]) 
    title(sprintf('r%0.3f, b%0.3f g%0.3f',rad2deg(atan(m))));
    
    figure;
    AnglePlot(subtheta(1),[0.8 0.1 0.1]) 
    AnglePlot(subtheta(2),[0.1 0.1 0.8]) 
    AnglePlot(subtheta(3),[0.1 0.8 0.1]) 
    title(sprintf('r%0.3f, b%0.3f g%0.3f',subtheta))

end
%add the angle between the lines and the rho value
rho = [];
rho(1) = corr(proj1(:,1),proj1(:,2));
rho(2) = corr(proj2(:,1),proj2(:,2));
rho(3) = corr(proj3(:,1),proj3(:,2));


end


function plotSubspaceAngles(xx)
fp = fig_params_cortdynamics;
% combined stats
x = cat(1,xx{:});
if size(x,2)>1
    theta = x(:,1)-x(:,2);
else
    theta = x; 
end

%get a bootstrapped distribution
[ci,thetaboot] = bootci(1000,@circ_mean, deg2rad(theta));
ci = rad2deg(ci);
thetaboot = rad2deg(thetaboot);
figure; hold on;
histogram(thetaboot,'facecolor',[0.35 0.35 0.35],'edgecolor','none','BinWidth',1)
yvals = get(gca,'ylim');
plot([0 0],yvals,'color',[0.25 0.25 0.25],'linewidth',1.5)
mu = rad2deg(circ_mean(deg2rad(thetaboot)));
err = rad2deg(circ_std(deg2rad(thetaboot),[],[],1));
plot([mu mu],yvals,'color',[0.8 0.1 0.1],'linewidth',1.5)
xval = get(gca,'xlim');
xlim([floor(xval(1)) max(ceil(xval(2)),5)]);
xlabel('\Delta theta (degrees)');
ylabel('bootstraps');
p = sum(thetaboot>=0)/numel(thetaboot);
title(sprintf('%0.2f +/- %0.2f \n CI: %0.2f %0.2f | p %0.3f',mu,err,ci(1),ci(2),p),'FontWeight','normal')
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 3 4],[10 10 14 10])

figure; hold on;
histogram(theta,'facecolor',[0.5 0.5 0.5],'edgecolor','none','BinWidth',5) %5 def
yvals = get(gca,'ylim');
plot([0 0],yvals,'color','k','linewidth',1)
mu = rad2deg(circ_mean(deg2rad(theta)));
err = rad2deg(circ_std(deg2rad(theta)))/sqrt(numel(theta));
plot([mu mu],yvals,'color',[0.8 0.1 0.1],'linewidth',1)
plot([ci(1) ci(1)],yvals,'color',[0.8 0.1 0.1],'linewidth',1,'linestyle','--')
plot([ci(2) ci(2)],yvals,'color',[0.8 0.1 0.1],'linewidth',1,'linestyle','--')
xlabel('\Delta theta (degrees)');
ylabel('count');
title(sprintf('%0.2f +/- %0.2f \n p %0.3f',mu,err,p),'FontWeight','normal')
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 3 4],[10 10 14 10])
set(gca,'ylim',yvals)


end


function x = GetSubspaceAngleForRecordings(data,cur_motif,area,cur_rec)

theta = cell(3,1);
thetacon = cell(3,1);

for i = 1:3
   [theta{i,1},thetacon{i,1}] = CompareToSubspace(data,cur_rec,cur_motif,[area(i),area(ismember(1:3,i)==0)],1000,i); 
end

theta = cellfun(@(x) deg2rad(x),theta,'UniformOutput',0);
thetacon = cellfun(@(x) deg2rad(x),thetacon,'UniformOutput',0);
%these numbers should be negative (i.e. smaller angle) 
x = cell(4,1);
%shared-->area 2 should be lower angle in B 
x(1,:) = cellfun(@(x,y) [circ_mean(y(2,:)),circ_mean(x(2,:))], theta(1,:),thetacon(1,:),'UniformOutput',0);
%shared-->area 1 should be lower angle in A
x(2,:) = cellfun(@(x,y) [circ_mean(x(1,:)),circ_mean(y(1,:))], theta(1,:),thetacon(1,:),'UniformOutput',0);
%area 1-->shared should be lower angle in A
x(3,:) = cellfun(@(x,y) [circ_mean(x(1,:)),circ_mean(y(1,:))], theta(2,:),thetacon(2,:),'UniformOutput',0);
%area 2-->shared should be lower angle in b
x(4,:) = cellfun(@(x,y) [circ_mean(y(1,:)),circ_mean(x(1,:))], theta(3,:),thetacon(3,:),'UniformOutput',0);
x = cellfun(@(x) rad2deg(x), x,'UniformOutput',0);

end


function x = GetSubspaceAngleAcrossRecordings(data,cur_motif,area,verbose)

if nargin <4; verbose =1; end
theta = cell(3,6);
thetacon = cell(3,6);
COUNT = 1;
for cur_rec = 1:6
    try
    for i = 1:3
        [theta{i,COUNT},thetacon{i,COUNT}] = CompareToSubspace(data,cur_rec,cur_motif,[area(i),area(ismember(1:3,i)==0)],1000,i); 
    end
    catch
    end
    COUNT = COUNT+1;
end
idx = cellfun(@(x) isempty(x), theta(1,:));
theta(:,idx)=[];
thetacon(:,idx)=[];

theta = cellfun(@(x) deg2rad(x),theta,'UniformOutput',0);
thetacon = cellfun(@(x) deg2rad(x),thetacon,'UniformOutput',0);
%these numbers should be negative (i.e. smaller angle) 
x = cell(4,sum(idx==0));
%shared-->area 2 should be lower angle in B 
x(1,:) = cellfun(@(x,y) [circ_mean(y(2,:)),circ_mean(x(2,:))], theta(1,:),thetacon(1,:),'UniformOutput',0);
%shared-->area 1 should be lower angle in A
x(2,:) = cellfun(@(x,y) [circ_mean(x(1,:)),circ_mean(y(1,:))], theta(1,:),thetacon(1,:),'UniformOutput',0);
%area 1-->shared should be lower angle in A
x(3,:) = cellfun(@(x,y) [circ_mean(x(1,:)),circ_mean(y(1,:))], theta(2,:),thetacon(2,:),'UniformOutput',0);
%area 2-->shared should be lower angle in b
x(4,:) = cellfun(@(x,y) [circ_mean(y(1,:)),circ_mean(x(1,:))], theta(3,:),thetacon(3,:),'UniformOutput',0);
x = cellfun(@(x) rad2deg(x), x,'UniformOutput',0);

if verbose ==1
rho = cell(3,6);
COUNT = 1;
for cur_rec = 1:6
    for i = 1:3
        try
            [~,~,y,yy,aa,bb,yyy] = RunPCA(data,cur_rec,cur_motif,area(i));
            rho{i,COUNT} = PlotExampleSubspaceRelationships(data,cur_rec,cur_motif,[area(i),area(ismember(1:3,i)==0)],y,yy,aa,bb,yyy,0); 
        catch
            rho{i,COUNT} = [NaN,NaN,NaN];
        end
    end
    COUNT = COUNT+1;
end

fp = fig_params_cortdynamics;
%get the difference in correlation between motif B to self and motif B to A
%subspace. 
r = cellfun(@(x) fisherZ(x),rho,'UniformOutput',0);
withinB = cellfun(@(x) abs(x(2)-x(3)),r,'UniformOutput',1);
ABself = cellfun(@(x) abs(x(1)-x(2)),r,'UniformOutput',1);
ABA = cellfun(@(x) abs(x(1)-x(3)),r,'UniformOutput',1);

% withinB = cellfun(@(x) (x(2)-x(3)),r,'UniformOutput',1);
% ABself = cellfun(@(x) (x(1)-x(2)),r,'UniformOutput',1);
% ABA = cellfun(@(x) (x(1)-x(3)),r,'UniformOutput',1);
r = [withinB(:),ABself(:),ABA(:)];
r(isnan(r(:,1)),:)=[];
[p,~,stats] = anova1(r,[],'off');
figure; hold on;       
boxplot(r,[],'notch','off')
f = get(get(gca,'children'),'children');
t = get(f,'tag');
idx = find(strcmp(t,'Box')==1);
for j = 1:numel(idx)
    set(f(idx(j)), 'Color', [0.25 0.25 0.25],'LineWidth',1.5);
end
idx = find(strcmp(t,'Median')==1);    
for j = 1:numel(idx)
    set(f(idx(j)), 'Color', 'k','LineWidth',1.5);
end
set(gca,'XTickLabel',{'Within','Between AB','Between AA'},'XTickLabelRotation',45)
fp.FormatAxes(gca);  box on; grid on; axis square; 
fp.FigureSizing(gcf,[3 3 2 3],[10 10 14 10])    
ylabel('\Delta Rho_z');
pp = multcompare(stats,'Display','off');
title(sprintf('%0.3f | %0.3f,%0.3f,%0.3f \n %0.2f+/-%0.2f \n %0.2f=/-%0.2f \n %0.2f=/-%0.2f',...
    p,pp(:,end),nanmean(r(:,1)),sem(r(:,1)),nanmean(r(:,2)),sem(r(:,2)),nanmean(r(:,3)),sem(r(:,3))),'fontweight','normal','fontsize',10);

%also plot a combined version with ttest
figure; hold on;  
grp = cat(1,ones(size(r,1),1),2*ones(2*size(r,1),1));
boxplot(r(:),grp,'notch','off')
f = get(get(gca,'children'),'children');
t = get(f,'tag');
idx = find(strcmp(t,'Box')==1);
for j = 1:numel(idx)
    set(f(idx(j)), 'Color', [0.25 0.25 0.25],'LineWidth',1.5);
end
idx = find(strcmp(t,'Median')==1);    
for j = 1:numel(idx)
    set(f(idx(j)), 'Color', 'k','LineWidth',1.5);
end
set(gca,'XTickLabel',{'Within','Between'},'XTickLabelRotation',45)
fp.FormatAxes(gca);  box on; grid on; axis square; 
fp.FigureSizing(gcf,[3 3 2 3],[10 10 14 10])    
ylabel('\Delta Rho_z');
p = ttest2(r(:,1),cat(1,r(:,2),r(:,2)),'tail','right');
title(sprintf('%0.3f \n %0.2f+/-%0.2f \n %0.2f=/-%0.2f',...
    p,nanmean(r(:,1)),sem(r(:,1)),nanmean(cat(1,r(:,2),r(:,3))),sem(cat(1,r(:,2),r(:,3)))),'fontweight','normal','fontsize',10);
end


end



function [theta_pca_in, theta_pca] = GetActivityAngleAcross(data,cur_motif,area,cur_rec)
a = cell(1,3);
b = cell(1,3);
theta = cell(1,3);
thetacon = cell(1,3);
for i = 1:3
    [~,~,y,yy] = RunPCA(data,cur_rec,cur_motif,area(i));
    a{i} = nanmean(y(:,1:3,:),3);
    b{i} = nanmean(yy(:,1:3,:),3);
end
% compare the angle between them
[theta_pca_in, theta_pca] = GetDiscrimitability(a,b); 
theta_pca_in = rad2deg(circ_mean(deg2rad(theta_pca_in)));
end

function CompareAngles(theta,thetacon,area)
figure;
[mu,err] = ViolinAnglePlot(theta{1}(2,:),[0.8 0.1 0.1],1);
[mu2,err2] = ViolinAnglePlot(thetacon{1}(2,:),[0.1 0.1 0.8],2);
set(gca,'xlim',[0.5 2.5])
set(gca,'xtick',[1,2]);
% figure;
% [mu,err] = AnglePlot(theta{1}(2,:),[0.8 0.1 0.1]);
% [mu2,err2] = AnglePlot(thetacon{1}(2,:),[0.1 0.1 0.8]);
temp = theta{1}(2,:)-thetacon{1}(2,:);
if sum(temp)>0; p = sum(temp<0)/numel(temp); else; p = sum(temp>0)/numel(temp); end
title(sprintf('Similarity between activity\nand subspace %s-->%s %0.3f\n Motif A (red) %s-->%s \nMotif B (blue) %s-->%s \n %0.2f +/- %0.2f | %0.2f +/- %0.2f',...
    area{1},area{3},p,area{1},area{2},area{1},area{3},mu,err,mu2,err2),'fontweight','normal','fontsize',8);

a = [1,2;2,1;3,1];
for i = 1:3
    figure;
    [mu,err] = ViolinAnglePlot(theta{i}(1,:),[0.8 0.1 0.1],1);
    [mu2,err2] = ViolinAnglePlot(thetacon{i}(1,:),[0.1 0.1 0.8],2);
    set(gca,'xlim',[0.5 2.5])
    temp = theta{i}(1,:)-thetacon{i}(1,:);
    if sum(temp)>0; p = sum(temp<0)/numel(temp); else; p = sum(temp>0)/numel(temp); end
    title(sprintf('Similarity between activity\nand subspace %s-->%s %0.3f\n Motif A (red) %s-->%s \nMotif B (blue) %s-->%s \n %0.2f +/- %0.2f | %0.2f +/- %0.2f',...
        area{a(i,1)},area{a(i,2)},p,area{1},area{2},area{1},area{3},mu,err,mu2,err2),'fontweight','normal','fontsize',8);
end

end

function [mu,err] = ViolinAnglePlot(theta,col,xpos)
fp = fig_params_cortdynamics;
CompareViolins(theta,fp,'plotspread',0,'divfactor',1,'xpos',xpos,'plotaverage',0,'col',col);
mu = rad2deg(circ_mean(deg2rad(theta)));
line([xpos-0.075 xpos+0.075],[mu mu],'Color',col,'LineWidth',3)
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 2.5 4],[10 10 14 10])
ylabel('Angle (degrees)');
err = rad2deg(circ_std(deg2rad(theta),[],[],2));
hold on; 
yval = get(gca,'ylim');
if yval(2)>90;
    set(gca,'ylim',[yval(1),90]);
end
end


function [mu,err] = AnglePlot(theta,col)    
fp = fig_params_cortdynamics;
edges = deg2rad(1:2:90);
theta = deg2rad(theta);
polarhistogram(theta,edges,'edgecolor','none','facecolor',col);    
set(gca,'ThetaLim',[0 90])
hold on;
rval = get(gca,'RLim');
mu = circ_mean(theta);
polarplot([mu, mu],rval,'linewidth',1.5,'color','k')
mu = rad2deg(mu);
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 4 4],[10 10 14 10])
set(gca,'GridAlpha',0.4)
err = rad2deg(circ_std(theta,[],[],2));
end


function [theta,thetacon,B,c] = CompareToSubspace(data,cur_rec,cur_motif,area,nboot,cur_i)
if nargin <5; nboot = 1000; end
rng('default');

%load the data
a = loadFunc(data,cur_rec,cur_motif(1),area{1},0); 
b = loadFunc(data,cur_rec,cur_motif(2),area{1},0); 

%relative to the subspaces of a and b
%if going x to y or y to x then use motif a subspaces
%CAMDEN - this looks like cur_i =1 and 2 are the same, but the area is
%changing with each call of this function so they are not. See above. 
if cur_i ==1
    B = loadCoef(data,cur_rec,cur_motif(1),area([1,2])); 
    BB = loadCoef(data,cur_rec,cur_motif(2),area([1,3]));   
elseif cur_i==2
    B = loadCoef(data,cur_rec,cur_motif(1),area([1,2])); 
    BB = loadCoef(data,cur_rec,cur_motif(2),area([1,3]));    
elseif cur_i==3
    B = loadCoef(data,cur_rec,cur_motif(2),area([1,2])); 
    BB = loadCoef(data,cur_rec,cur_motif(1),area([1,3]));
end

theta = NaN(1,nboot); %each row is a comparison with a subspace
x = reshape(a,[size(a,1),size(a,2)*size(a,3)])';
x = x-nanmean(x); %mean center
coef = pca(x);
c = coef(:,1); %save off for later
%angle between local representation and the subspace representation  
if nboot >1
    idx = arrayfun(@(n) datasample(1:size(a,1),size(a,1)),1:nboot,'UniformOutput',0);
else
    idx  = {1:size(a,1)};
end
for j = 1:nboot %bootstrap the theta distribution
    theta(1,j) = rad2deg(subspace(coef(idx{j},1),B(idx{j},1)));
    theta(2,j) = rad2deg(subspace(coef(idx{j},1),BB(idx{j},1)));
end

thetacon = NaN(1,nboot);
x = reshape(b,[size(b,1),size(b,2)*size(b,3)])';
x = x-nanmean(x); %mean center
coef = pca(x);
%angle between local representation and the subspace representation  
for j = 1:nboot %bootstrap     
    thetacon(1,j) = rad2deg(subspace(coef(idx{j},1),B(idx{j},1)));
    thetacon(2,j) = rad2deg(subspace(coef(idx{j},1),BB(idx{j},1)));
end
c = [c,coef(:,1)];
end


function [theta_in,theta_join] = GetDiscrimitability(a,b)
    n = numel(a);
    theta_in = NaN(1,n);
    for i = 1:n
        theta_in(i) = rad2deg(subspace(a{i},b{i}));
    end        
    x = cat(3,a{:});
    x = squeeze(x(:,1,:));
    y = cat(3,b{:});
    y = squeeze(y(:,1,:));
    theta_join = rad2deg(subspace(x,y));
end


function B = loadCoef(data,cur_rec,cur_motif,area_name) 
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    paired_areas = data{cur_rec}(cur_motif).paired_areas;
    
    a = find(strcmp(area_label,area_name{1}));
    b = find(strcmp(area_label,area_name{2}));
    idx = find(ismember(paired_areas,[a,b],'rows')); 
    
    B = data{cur_rec}(cur_motif).rrr_B{idx};
end



function [s,sig] = CompareActivity(data,cur_rec,cur_motif,area,verbose)
if nargin <5; verbose = 1; end
fp = fig_params_cortdynamics;
%get all motifs
rng('default');
col = [0.1 0.1 0.8; 0.8 0.1 0.1];
for i = 1:numel(area)
    %compare the values to subtracted baseline (yes, this is what we
    %shoudl've done from the start - easier than the division baseline
    abase = StatCompareToBaselin(data,cur_rec,cur_motif(1),area{i}); 
    bbase = StatCompareToBaselin(data,cur_rec,cur_motif(1),area{i});
    abase = nanmean(squeeze(max(abase,[],2)),2);
    bbase = nanmean(squeeze(max(bbase,[],2)),2);
    [~,pabase] = ttest(abase,0,'tail','right');
    [~,pbbase] = ttest(bbase,0,'tail','right');
    
    a = loadFunc(data,cur_rec,cur_motif(1),area{i},0); 
    b = loadFunc(data,cur_rec,cur_motif(2),area{i},0); 
    %get average peak activity across neurons per area
    a = nanmean(squeeze(max(a,[],2)),2);
    b = nanmean(squeeze(max(b,[],2)),2);
 
    if verbose
    %simple mean + sem plot
    figure; hold on; 
    plot([1,2],[nanmean(a),nanmean(b)],'color',[0.5 0.5 0.5],'linewidth',1)
    errorbar(1,nanmean(a),sem(a),'marker','o','linewidth',1,'color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',fp.markersizesmall);
    errorbar(2,nanmean(b),sem(b),'marker','o','linewidth',1,'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'MarkerSize',fp.markersizesmall);
    
    [~,p] = ttest(a,b);
    title(sprintf('M %d & %d %s \n rec %d p=%0.4f \n %0.3f+/-%0.3f\n %0.3f+/-%0.3f\n %0.3f %0.3f baseline',...
        cur_motif(1),cur_motif(2),area{i},cur_rec,p,nanmean(a),sem(a),nanmean(b),nanmean(b),pabase,pbbase),'fontweight','normal');    
    set(gca,'XTickLabel',cur_motif,'xlim',[0.5 2.5])
    yval = get(gca,'ylim');
    set(gca,'ylim',[round(yval(1)-0.025,2),round(yval(2)+0.025,2)]);
    fp.FormatAxes(gca);  box on; grid on; axis square; 
    fp.FigureSizing(gcf,[3 3 1.25 2],[10 10 14 10])   
    yval = get(gca,'ylim');
    set(gca,'ytick',yval);
    ylabel('Activity');    
    %boxplot  
    grp = cat(1,ones(numel(a),1),ones(numel(b),1)*2);
    figure; hold on        
    boxplot(cat(1,a,b),grp,'notch','on')
    f = get(get(gca,'children'),'children');
    t = get(f,'tag');
    idx = find(strcmp(t,'Box')==1);
    for j = 1:numel(idx)
        set(f(idx(j)), 'Color', col(j,:),'LineWidth',1.5);
    end
    idx = find(strcmp(t,'Median')==1);    
    for j = 1:numel(idx)
        set(f(idx(j)), 'Color', 'k','LineWidth',1.5);
    end
    set(gca,'XTickLabel',cur_motif)
    fp.FormatAxes(gca);  box on; grid on; axis square; 
    fp.FigureSizing(gcf,[3 3 3 3],[10 10 14 10])    
    ylabel('Activity');
    
    idx = a>b; 
    plot([1.25,1.75],[a(idx),b(idx)],'color',[0.8 0.1 0.1],'linewidth',.5,'linestyle','--')
    plot([1.25,1.75],[a(idx==0),b(idx==0)],'color',[0.1 0.1 0.8],'linewidth',1,'linestyle',':')        
    ylim([0 ceil(max(cat(1,a(:),b(:))))])
    p = signrank(a,b);
    title(sprintf('M %d & %d %s \n rec %d p=%0.2d',cur_motif(1),cur_motif(2),area{i},cur_rec,p),'fontweight','normal');    
    end
    if nanmean(a)>nanmean(b)
        s(i)=1;        
    else
        s(i)=0;       
    end
    sig(i) = ttest2(a,b);

end
s = s(2:end);
sig = sig(2:end);
    
end

function PlotJointProjections(a,b,cur_motif,area,cur_rec)
fp = fig_params_cortdynamics;

for i = 1
   %just the patch
   x = cat(3,a{:});
   y = cat(3,b{:});
   x = squeeze(x(:,i,:));   
   y = squeeze(y(:,i,:));
   figure; hold on; 
   x = movmean(x,2,1);
   y = movmean(y,2,1);
   [X,Y,Z] = Get3DPlane(x,0.75);
   patch(X,Y,Z,[0.8 0.1 0.1],'EdgeColor','none','FaceAlpha',.5);
   [X,Y,Z] = Get3DPlane(y,0.75);
   patch(X,Y,Z,[0.1 0.1 0.8],'EdgeColor','none','FaceAlpha',.5);   
   xlabel(sprintf('PC %d | %s',i,area{1}))
   ylabel(sprintf('PC %d | %s',i,area{2}))
   zlabel(sprintf('PC %d | %s',i,area{3}))
   set(gca,'CameraPosition',[9 3 1.7])
   fp.FormatAxes(gca);  box on; grid on; axis square; 
   fp.FigureSizing(gcf,[3 3 3 3],[10 10 14 10])
   title(sprintf('M %d & %d |rec %d',cur_motif(1),cur_motif(2),cur_rec),'fontweight','normal'); 
   
   
   x = cat(3,a{:});
   y = cat(3,b{:});
   x = squeeze(x(:,i,:));   
   y = squeeze(y(:,i,:));
   figure; hold on; 
   x = movmean(x,2,1);
   y = movmean(y,2,1);
   [X,Y,Z] = Get3DPlane(x,0.75);
   patch(X,Y,Z,[0.8 0.1 0.1],'EdgeColor','none','FaceAlpha',.2);
   [X,Y,Z] = Get3DPlane(y,0.75);
   patch(X,Y,Z,[0.1 0.1 0.8],'EdgeColor','none','FaceAlpha',.2);   
   plot3(x(:,1),x(:,2),x(:,3),'marker','o','linestyle','-','linewidth',1,'color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',fp.markersizesmall)
   plot3(y(:,1),y(:,2),y(:,3),'marker','o','linestyle','-','linewidth',1,'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'MarkerSize',fp.markersizesmall)
   
   %plot a marker at the start
   plot3(x(1,1),x(1,2),x(1,3),'marker','d','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)
   plot3(y(1,1),y(1,2),y(1,3),'marker','d','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)
   xlabel(sprintf('PC %d | %s',i,area{1}))
   ylabel(sprintf('PC %d | %s',i,area{2}))
   zlabel(sprintf('PC %d | %s',i,area{3}))
   set(gca,'CameraPosition',[9 3 1.7])
   fp.FormatAxes(gca);  box on; grid on; axis square; 
   fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
   title(sprintf('M %d & %d |rec %d',cur_motif(1),cur_motif(2),cur_rec),'fontweight','normal');
   
   %no patch version
   x = cat(3,a{:});
   y = cat(3,b{:});
   x = squeeze(x(:,i,:));   
   y = squeeze(y(:,i,:));
   figure; hold on; 
   x = movmean(x,2,1);
   y = movmean(y,2,1); 
   plot3(x(:,1),x(:,2),x(:,3),'marker','o','linestyle','-','linewidth',1,'color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',fp.markersizesmall)
   plot3(y(:,1),y(:,2),y(:,3),'marker','o','linestyle','-','linewidth',1,'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'MarkerSize',fp.markersizesmall)
   
   %plot a marker at the start
   plot3(x(1,1),x(1,2),x(1,3),'marker','d','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)
   plot3(y(1,1),y(1,2),y(1,3),'marker','d','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)
   xlabel(sprintf('PC %d | %s',i,area{1}))
   ylabel(sprintf('PC %d | %s',i,area{2}))
   zlabel(sprintf('PC %d | %s',i,area{3}))
   set(gca,'CameraPosition',[9 3 1.7])
   fp.FormatAxes(gca);  box on; grid on; axis square; 
   fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
   title(sprintf('M %d & %d |rec %d',cur_motif(1),cur_motif(2),cur_rec),'fontweight','normal');   
end


end



function PlotProjection(y,yy,cur_motif,area_name,cur_rec,comp)
fp = fig_params_cortdynamics;

if nargin <6; comp = [1,2;1,3;2,3]; end
for i =1:size(comp,1)
    a = nanmean(y,3);
    b = nanmean(yy,3);    
    a = [a(:,comp(i,1)),a(:,comp(i,2))];
    b = [b(:,comp(i,1)),b(:,comp(i,2))];
    figure; hold on; 
    plot(a(:,1),a(:,2),'marker','o','linestyle','-','linewidth',1,'color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',fp.markersizesmall)
    plot(b(:,1),b(:,2),'marker','o','linestyle','-','linewidth',1,'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'MarkerSize',fp.markersizesmall)
    %plot a marker at the start
    plot(a(1,1),a(1,2),'marker','o','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)
    plot(b(1,1),b(1,2),'marker','o','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)    
    xval = cat(1,a(:,1),b(:,1));
    AddLSline(a(:,1),a(:,2),xval,[0.8 0.1 0.1],1.5);
    AddLSline(b(:,1),b(:,2),xval,[0.1 0.1 0.8],1.5);
    fp.FormatAxes(gca);  box on; grid on; axis square; 
    if i > 1
        fp.FigureSizing(gcf,[3 3 2 2],[10 10 10 10])
    else
        fp.FigureSizing(gcf,[3 3 4 4],[10 10 10 10])
    end
    xlabel(sprintf('PC %d',comp(i,1)));
    ylabel(sprintf('PC %d',comp(i,2)));
    title(sprintf('M %d & %d | %s |rec %d',cur_motif(1),cur_motif(2),area_name{1},cur_rec),'fontweight','normal');
end



end




function x = StatCompareToBaselin(data,cur_rec,cur_motif,area_name)
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    x = area_val{strcmp(area_label,area_name)};

    %to compare with baseline activity
    x = normalizeToBaseline(x,[1:2],'meansubtract');

end


function x = loadFunc(data,cur_rec,cur_motif,area_name,flag)
    if nargin <5; flag = 1; end
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    x = area_val{strcmp(area_label,area_name)};

    %normalize to baseline
%     x = normalizeToBaseline(x,[1:2],'mean');
    
    %to compare with baseline activity
    x = normalizeToBaseline(x,[1:2],'meansubtract');

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


function [coef,coef2,y,yy,a,b,yyy] = RunPCA(data,cur_rec,cur_motif,area)
% if nargin <5; shufflag = 1; end
rng('default');

%get all motifs
a = loadFunc(data,cur_rec,cur_motif(1),area{1},0); 
b = loadFunc(data,cur_rec,cur_motif(2),area{1},0); 

%% uncomment to perform on withheld
%get random subset of data from each motif
% n = min(size(a,3),size(b,3));
% idx = randperm(size(a,3),n); %so same random trials in x and y
% a = a(:,:,idx);
% idx = randperm(size(b,3),n); %so same random trials in x and y
% b = b(:,:,idx);
% testidx = randperm(n,floor(0.4*n));
% trainidx = ismember(1:n, testidx)==0;

%% fit the PCA on a each motif
% x = a(:,:,trainidx);
x = a;
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
x = x-nanmean(x); %mean center
coef = pca(x);

% x = b(:,:,trainidx);
x = b; 
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
x = x-nanmean(x); %mean center
coef2 = pca(x);


% a = a(:,:,testidx);
% b = b(:,:,testidx);

%% projection into the space of one of the motifs
y = cell(numel(size(a,3)),1);
for i = 1:size(a,3)
    x = a(:,:,i)';
    x = x-nanmean(x);
    y{i} = x*coef(:,1:4);
end
yy = cell(numel(size(b,3)),1); %B on it's motif 1s
for i = 1:size(b,3)
    x = b(:,:,i)';
    x = x-nanmean(x);
    yy{i} = x*coef(:,1:4);    
end
yyy = cell(numel(size(b,3)),1); %B on it's own PC
for i = 1:size(b,3)
    x = b(:,:,i)';
    x = x-nanmean(x);
    yyy{i} = x*coef2(:,1:4);    
end
y = cat(3,y{:});
yy = cat(3,yy{:});
yyy = cat(3,yyy{:});



end


function [theta,thetacon] = RunPCAacrossTime(data,cur_rec,cur_motif,area)
% if nargin <5; shufflag = 1; end
warning('needs to be changed to circular statistics');
rng('default');

%load the data
a = loadFunc(data,cur_rec,cur_motif(1),area{1},0); 

%relative to the subspaces of a and b
B = loadCoef(data,cur_rec,cur_motif(1),area([1,2])); 
BB = loadCoef(data,cur_rec,cur_motif(2),area([1,3])); 

theta = NaN(size(a,2),1000);
thetacon = NaN(size(a,2),1000);
for i = 1:size(a,2)
    x = a(:,i,:);
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    x = x-nanmean(x); %mean center
    coef = pca(x);
    %angle between local representation and the subspace representation     
    for j = 1:1000 %bootstrap
        idx = datasample(1:size(coef,1),size(coef,1));
        theta(i,j) = rad2deg(subspace(coef(idx,1),B(idx,1)));
        thetacon(i,j) = rad2deg(subspace(coef(idx,1),BB(idx,1)));
    end
end

figure; hold on; 
ci = cat(1,prctile(theta',2.5),prctile(theta',97.5))';
ci(:,1) = nanmean(theta,2)-ci(:,1);
ci(:,2) = ci(:,2)-nanmean(theta,2);
shadedErrorBar(1:12,nanmean(theta,2),ci,'lineprops',{'color',[0.8 0.1 0.1 0.75],'linewidth',2});
ci = cat(1,prctile(thetacon',2.5),prctile(thetacon',97.5))';
ci(:,1) = nanmean(thetacon,2)-ci(:,1);
ci(:,2) = ci(:,2)-nanmean(thetacon,2);
shadedErrorBar(1:12,nanmean(thetacon,2),ci,'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});




end






% 
% 
% n = size(a{1},1);
% predictors = cat(1,a{1},b{1});
% response = cat(1,ones(n,1),zeros(n,1));
% %fit a hyperplane
% classificationSVM = fitcsvm(predictors,response,'KernelFunction','linear');
% 
% %classification plane in 3D
% x = cellfun(@(x) x(:,1), a,'UniformOutput',0);
% y = cellfun(@(x) x(:,1), b,'UniformOutput',0);
% x = cat(2,x{:});
% y = cat(2,y{:});
% 
% predictors = cat(1,x,y);
% response = cat(1,ones(n,1),2*ones(n,1));
% 
% nfold = n*2;
% c = cvpartition(n*2,'KFold',nfold);  %loocv 
% predictors = cat(1,a{1},b{1});
% response = cat(1,ones(n,1),2*ones(n,1));
% rng('default')
% pred_response = NaN(nfold,1);
% for i = 1:nfold
%     classificationSVM = fitcsvm(predictors(c.training(i),:),response(c.training(i),:),'KernelFunction','linear');
%     %predict test data
%     svmPredictFcn = @(x) predict(classificationSVM, x);
%     validationPredictFcn = @(x) svmPredictFcn(x);
%     [pred_response(i), ~] = validationPredictFcn(predictors(c.training(i),:));        
% end
% true = arrayfun(@(n) response(c.test(n)),1:nfold);
% pred_response==true';


% function theta_perm = ShuffledTheta(a,b)
% rng('default');
% nperm =1000; 
% theta_perm = NaN(1,nperm); 
% for cur_p = 1:nperm
%     aa = a;
%     bb = b; 
%     for j = 1:3 %break independences of a and b by mixing
%         temp = cat(1,a{j},b{j});
%         temp = temp(:);
%         idx = randperm(size(temp,1),size(temp,1)/2);
%         aa{j} = reshape(temp(idx),12,3);
%         temp(idx)=[];
%         bb{j} = reshape(temp,12,3);
%     end
%     [~, theta_perm(cur_p)] = GetDiscrimitability(aa,bb);
% end
% end
% 
% 
% function theta_perm = TrialShuffledTheta(data,cur_rec,cur_motif,area)
% rng('default');
% nperm =100; 
% theta_perm = NaN(1,nperm); 
% for cur_p = 1:nperm
%     a = cell(1,3);
%     b = cell(1,3);
%     for i = 1:3
%         [~,~,y,yy] = RunPCA(data,cur_rec,cur_motif,area(i),cur_p*i);
% %         PlotProjection(y,yy,cur_motif,area(i),cur_rec,[1,2]);
%         a{i} = nanmean(y(:,1:3,:),3);
%         b{i} = nanmean(yy(:,1:3,:),3);
%     end
%     [~, theta_perm(cur_p)] = GetDiscrimitability(a,b);
% end
% end

% %permutation test to break the difference between a and b by merging and randomly sampling trials
% if shufflag>1
%     rng(shufflag)
%     temp = cat(3,a,b);
%     idx = randperm(size(temp,3),size(a,3));
%     a = temp(:,:,idx);
%     temp(:,:,idx)=[];
%     b = temp; 
% end











