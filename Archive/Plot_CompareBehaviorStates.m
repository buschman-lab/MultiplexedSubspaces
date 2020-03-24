function fh = Plot_CompareBehaviorStates(weight,num_states,yvals)

fp = fig_params;
sig_color = [0.75 0.75 0.75; 0.25 0.25 0.25];
col = getColorPalet(num_states);
pval = [];
fstat = [];
y=[];
sig_comparison = {};
for i = 1:size(weight,1)
   temp = cat(1,weight{i,:})';
   temp = (temp-nanmean(temp(:)))/nanstd(temp(:));
   y(i,:) = nanmean(temp);
   [pval(i),tab,stats] = anova1(temp,[],'off');
   c = multcompare(stats,'Display','off');
   temp_tukey =c(:,1:2);
   sig_comparison{i} = temp_tukey(sum(c(:,end)<0.05,2)>0,:);   
   fstat(i) = tab{2,5};
end

figure('Position',[0 0 1000 1000]); hold on; 
s1 = subplot(312,'Units','centimeters','Position',[8 6 10 4.5]); hold on
s2 = subplot(311,'Units','centimeters','Position',[8 11 10 1.5]); hold on

axes(s1);
y = y';
imagesc(y,[yvals])
colormap(gca,flipud(redgreencmap(64,'Interpolation','linear')));
c = colorbar;
set(c,'YTick',[yvals(1), 0, yvals(2)]);
ylabel(c,{'Motif Intensity';'(zscore)'},'FontSize',16,'Fontweight','normal','FontName','Arial');
set(c,'units','centimeters','position',[18.25 6 0.5 4.5])

for x_grid = 0.5:1:size(y,2)+0.5
    line([x_grid,x_grid],[0.5,size(y,1)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid = 0.5:1:size(y,1)+0.5
    line([0.5,size(y,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  
xlabel('Basis Motifs')


xlim([0.5 size(y,2)+.5])
ylim([0.5 size(y,1)+0.5])
set(gca,'YColor','none')
setFigureDefaults
% For each significant comparison put a box around motifs that were different
for i = 1:numel(sig_comparison)
   temp = unique(sig_comparison{i}(:));
   for j = 1:numel(temp)
       plot(i,temp(j),'.','markersize',25,'color',[1,0.8,0.6])     
   end    
end

axes(s2); hold on
%Plot the fscore and the significance
plot(fstat,'LineWidth',2,'Marker','.','MarkerSize',5,'MarkerEdgeColor',[0.4 0.4 0.4],'Color',[0.4 0.4 0.4])
for i = 1:numel(pval)
   AddSig(1,pval(i),[i-0.1,i-0.1,fstat(i),fstat(i)],1,5,1,90)
end
%Change marker color for significant motifs (lighter)
is_sig = cellfun(@(x) ~isempty(x),sig_comparison,'UniformOutput',0);
scatter((1:1:size(weight,1)),fstat,50,sig_color(cell2mat(is_sig)+1,:),'filled')
xlim([0.5 size(weight,1)+.5])
set(gca,'XTick',(1:2:size(weight,2)),'Ytick',[min(get(gca,'Ytick')),max(get(gca,'Ytick'))]);
ylabel({'F';'stat';''},'Rotation',0,'Units','Centimeters','position',[11.25 1.4]);
set(gca,'yaxislocation','right')
set(gca,'XColor','none');
setFigureDefaults


set(gcf,'Position',[680   150  875  650]);

fh = gcf;
end