function Figure2f(behav_rho)
fp = fig_params_cortdynamics;

% Plot
rho = fisherZ(behav_rho); 
figure; hold on; 
[a,~] = pairedBootstrap(rho,@nanmean);
CompareViolins(a',fp,'plotspread',0,'divfactor',1,'plotaverage',1,'col',repmat({[0.15 0.15 0.15]},1,10),'distWidth',0.75,'connectline',[0.75 0.75 0.75]);
xlim([0.5 10.5])
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 5 4],[2 10 10 10])
xlabel('Subspace dimension');
ylabel('rho_z')
ylim([0.05 0.35])

end