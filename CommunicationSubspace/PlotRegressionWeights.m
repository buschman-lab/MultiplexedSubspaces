function fh = PlotRegressionWeights(b,m,reg,area_lbl,fp)
if nargin <5; fp = fig_params_cortdynamics; end

%full beta weights
figure; hold on;  
x = squeeze(b(1,m,:));
x(isnan(x))=[];
bar(x,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5,'EdgeAlpha',0);
lbl = squeeze(reg(1,m,:));
lblu = unique(lbl(~isnan(lbl)));
c = getColorPalet(numel(lblu));
p = cell(1,numel(lblu));
for i = 1:numel(lblu)
   p{i} = plot(find(lbl==lblu(i)),zeros(1,sum(lbl==lblu(i))),'color',c(i,:),'linewidth',2);
end
legend([p{:}],area_lbl(lblu),'location','bestoutside');
xlabel('nueron')

fp.FormatAxes(gca); 
box on;
fp.FigureSizing(gcf,[3 2 16 4],[10 10 26 8])
xlim([-10,size(x,1)+10])
fh = get(gcf);

end %function end