function PlotExampleEphys(cur_rec,tp)
fp = fig_params_cortdynamics;
%camden macdowell
if nargin <3; tp = [100,100+(60*15*1)]; end
[area_val,area_label] = GetExampleSpikingData(cur_rec,1);

%plot by area so you can mix and rearrange later
x = cat(1,area_val{:});
x = x(:,[tp(1):tp(2)]);
figure; hold on; 
imagesc(x,[0 1]); colorbar
set(gcf, 'Renderer', 'painters')
colormap(gca,flipud(gray));
xlim([-5,size(x,2)+0.5])
ylim([0,size(x,1)+0.5])

for i = 1:numel(area_label)    
    if i>1
        temp = cat(1,area_val{1:i-1});
        idx = [size(temp,1)+1,size(temp,1)+size(area_val{i},1)];
    else
        idx = [1,size(area_val{i},1)];
    end
    plot([-5 -5],idx,'color',fp.c_area(i,:),'linewidth',3); 
    text(-10, idx(1)+(idx(2)-idx(1))/2, area_label{i},'FontWeight','bold','HorizontalAlignment','right','Color',fp.c_area(i,:),'fontsize',fp.font_size,'FontName',fp.font_name)
    
end
set(gca,'YAxisLocation','right');
set(gca,'xtick',0:(15*10):(tp(2)-tp(1)),'xticklabel',(0:(15*10):(tp(2)-tp(1)))/15)
xlabel('time (seconds)');
ylabel('neurons')
fp.FormatAxes(gca); box on
set(gca,'Clipping','off')


%detailed





end %function end
