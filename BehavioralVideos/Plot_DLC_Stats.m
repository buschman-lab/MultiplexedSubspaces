function Plot_DLC_Stats(in_fns)

if nargin <1
   in_fns = GrabFiles('.csv');
end

[~,~,raw_data] = xlsread(in_fns{1});

raw_data = cell2mat(raw_data);

figure('position',[1105 524 433 454]); hold on; 
plot(raw_data(:,1),raw_data(:,2),'Linewidth',2,'Color',[0.5 0.5 0.5]);
title({'Markerless Tracking Was';'Trained Until Loss Plateaued'},'FontName','Arial','Fontweight','normal','FontSize',16); 
xlabel('Training Iteration');
ylabel('Loss');
ax = gca; 
ax.XAxis.Exponent = 0;
setFigureDefaults
set(gca,'position',[3 3 6 6])


handles = get(groot, 'Children');
saveCurFigs(handles,'-svg','DLC_Metrics',pwd,1);
close all;


