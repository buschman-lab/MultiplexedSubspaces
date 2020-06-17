function fh = DetectDictalEvents(fn)
%As in mussal et al., 2020 and steinmetz et al., 2017

temp = load(fn,'dff','opts');
signal = nanmean(conditionDffMat(temp.dff),2);
[~,~,widths,prominences] = findpeaks(signal,opts.fps); 

[~, name] = fileparts(fn);
mouse = name(1:regexp(name,'_','once')-1); %assumes Mouse#_ structure

figure; plot(widths,prominences,'.','markersize',0.1,'color','k');
xlim([0 0.75])
ylim([0 4]);
title(['Mouse ', mouse],'FontSize',16,'Fontweight','normal','FontName','Arial');
xlabel('peak width (s)');
ylabel('prominence (df/f)');
setFigureDefaults();

fh = gcf;

end


