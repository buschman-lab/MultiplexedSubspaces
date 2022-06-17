function Plot_SubspaceMultiplexing(data)
%Camden MacDowell - timeless
%looks trial-to-trial on the relationship between subspace dimensions
fp = fig_params_cortdynamics;

%get example data all recordings
rsq = cell(1,6);
for i = 1:6
   [~,~,~,~,rsq{i}] = SubspaceDim_Trial(data{i},0);
end
[~,area_label] = LoadVariable(data,'rel_performance',[]);

%adjust for missing
temp = cell(8,14);
temp([1,2,3,5,7,8],:) = rsq{3};
rsq{3} = temp;
temp = cell(8,14);
temp([1,2,3,5,6,7,8],:) = rsq{4};
rsq{4}=temp;

%% plot the average activity 
close all;
for a = 1:8
    x = rsq{1}(a,6); %one area across all trials for one motif on a given recording
    x = cat(1,x{:});
    x = sort(x,2,'descend');
    x =  cumsum(x,2)./sum(x,2);
    col = fp.c_area; col = col(strcmp(area_label,area_label{a}),:);
    figure; hold on;
    arrayfun(@(n) plot(1:size(x,2), x(n,:),'linestyle','-','marker','none','color',[0.85 0.85 0.85],'linewidth',0.5,'markersize',4), 1:size(x,1));
    plot(1:size(x,2), nanmean(x),'linestyle','-','marker','o','color',col,'linewidth',1.5,'markersize',5,'MarkerFaceColor','none');  
    % shadedErrorBar(1:size(x,2),nanmean(x),std(x,1),'lineprops',{'color',[col,0.25],'linewidth',2});
    xlim([1 size(x,2)]); ylim([0 max(x(:))])
    plot([0 size(x,2)],[0.8 0.8] ,'linestyle','--','color','k','linewidth',1.5);
    xlabel('# of dimensions (sorted)');
    ylabel({'Contribution to trial (r^2)'});
    title(sprintf('rec 1\n %s | motif %d',area_label{a},5),'fontweight','normal')
    fp.FormatAxes(gca);   box on; grid on
    fp.FigureSizing(gcf,[3 2 3.5 3.5],[10 10 10 10])
    set(gca,'xtick',[1,5,10,15]);
end

% After the example, plot for each areas
figure; hold on;
d = NaN(1,8);
for a = 1:8   
    x = cat(3,rsq{:});
    x = squeeze(x(a,:,:)); %one area across motifs
    x = cat(1,x{:});
    x = sort(x,2,'descend');
    x =  cumsum(x,2)./sum(x,2);
    d(a) = find(nanmean(x,1)>0.8,1,'first');
    col = fp.c_area; col = col(strcmp(area_label,area_label{a}),:);
    shadedErrorBar(1:size(x,2),nanmean(x,1),sem(x,1),'lineprops',{'color',[col 0.75],'linewidth',2});
end
set(gca,'xtick',[1,5,10,15]);
xlim([1 size(x,2)]); ylim([0 max(x(:))])
plot([0 size(x,2)],[0.8 0.8] ,'linestyle','--','color','k','linewidth',1.5);
xlabel('# of dimensions (sorted)');
ylabel({'Contribution to trial (r^2)'});
title(sprintf('%s | %d-%d',area_label{a},min(d),max(d)),'fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 3.5 3.5],[10 10 10 10])


%% projection of each trial along the first and second subspace




%% same for third and fourth





%% plot the PEV of SS1 for SS2 and SS3 versus SS4
all_areas = data{1}.area_label;
yy = cell(1,6);
for i = 1:6
   yy{i} = loadProjections(data{3},all_areas);
end
y = cat(6,yy{:});


%%
figure; hold on;
for a = 1:8  
    x = squeeze(y(a,:,:,:,:,:)); %one area across motifs/recs/trials
    x = permute(x,[3,2,4,1,5]);
    %just plot an example motif/rec
    x = (x(:,:,:,5,1));    
    figure; hold on;
    x = x(:,:,squeeze(~isnan(x(1,1,:))));
    z = nanmean(squeeze(x(1,:,:)),1);
    zz = nanmean(squeeze(x(2,:,:)),1);
    z = (z(:));
    zz = (zz(:));    
        
    plot(z,zz,'.','color',col,'markersize',2,'linestyle','none')
    AddLSline(z,zz,z,col);
    xlabel('projection along dim 1');
    ylabel('projection along dim 2');
    
end
    
%     x = squeeze(y(a,:,:,:,:,:)); %one area across motifs/recs/trials
%     x = permute(x,[2,3,4,1,5]);
%     %just plot an example motif/rec
%     x = x(:,:,:,5,6);
%    
%     col = fp.c_area; col = col(strcmp(area_label,area_label{a}),:);    
%     x = x(:,:,squeeze(~isnan(x(1,1,:))));
%     m = NaN(1,size(x,3));
%     for i = 1:size(x,3)
%         %get the slope
%         lm=fitlm(x(:,1,i)-mean(x(:,3,i)),x(:,4,i)-mean(x(:,1,i)));             % a linear model
%         mm=lm.Coefficients.Estimate;
%         m(i) = mm(2);        
%     end
%     m = abs(m);  %correct for reversal
%     %plot all the lines
%     figure; hold on; 
%     for i = 1:numel(m)
%         f = @(x) m(i)*x;
%         plot(f([0 1]),'linewidth',0.2,'color',[0.25 0.25 0.25])
%     end
%     
% end

plot(x(:,1),x(:,2),'.','color',col,'markersize',1,'linestyle','none')
end %function

%%


function yy = loadProjections(data,all_area)
ndim = 4;
rng('default')
area_label = data(1).area_label;
yy = NaN(8,14,12,ndim,1000);
% yy = NaN(8,14,ndim,1000);
%loop through each area
for cur_area = 1:numel(area_label) 
    fprintf('\n\t working on area %d of %d',cur_area,numel(area_label));       
    for cur_motif = 1:size(data,2)
        x = cat(1,data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==0});

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end,:);

        %remove the psth
        x = x-nanmean(x,3);

        %full model
        B = data(cur_motif).rrr_B{cur_area};        
        for i = 1:ndim
            if sum(B(:,i)>0)<sum(B(:,i)<0)
                B(:,i) = B(:,i)*-1;
            end  
        end
%         V = data(cur_motif).rrr_V{cur_area};
              
        %projection along each dimension
        for t = 1:size(x,3)
           xx = squeeze(x(:,:,t))'; 
           yy(strcmp(all_area,area_label{cur_area}),cur_motif,:,1:4,t) = xx*B(:,1:ndim);
%            yy(strcmp(all_area,area_label{cur_area}),cur_motif,1:4,t) = nanmean(xx*B(:,1:ndim),1);
%            for cur_d = 1:ndim
%               yy(strcmp(all_area,area_label{cur_area}),cur_motif,cur_d,t) = nanmean(xx*B(:,cur_d)*V(:,cur_d)',[1,2]);
%            end
        end  
        
    end    
end

%average within trials
% yy = squeeze(nanmean(yy,3));


end %function end





%% plot stability across recordings and areas
% rho = SubspaceDim_Trial(data{1},0);
% rho = fisherZ(rho);
% x = pairedBootstrap(rho',@nanmean);
% [~,ind] = sort(nanmean(x),'descend');
% figure; hold on;       
% col = fp.c_area;
% col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
% %sort
% col = col(ind);
% CompareViolins(x(:,ind)',fp,'plotspread',0,'divfactor',0.5,'plotaverage',1,'col',col,'distWidth',0.75);
% set(gca,'XTickLabel',area_label(ind),'XTickLabelRotation',45)
% yval = get(gca,'ylim');
% ylabel('% Generalization');
% fp.FormatAxes(gca);  box on; grid on
% fp.FigureSizing(gcf,[3 3 6 3],[10 10 10 10])
%% get the bootstrapped stats for consistency of this pattern across recordings/motifs





