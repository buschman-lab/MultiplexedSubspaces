function OrderActivityBySubspaceDimension(data,cur_rec)
%Camden MacDowell - timeless
tp = [1200 1200+(60*15*2)-1];
[area_val,area_label] = GetExampleSpikingData(cur_rec,1);

%load an example subspace (i.e. vis and motif 5)
area = 'VIS';
idx = strcmp(area_label,area);
motif = 5;
B = data{cur_rec}(motif).rrr_B{strcmp(area_label,area)};
V = data{cur_rec}(motif).rrr_V{strcmp(area_label,area)};

%get a 5 minute chunk of data
x = cat(1,area_val{idx==0});
y = area_val{idx==1};
x = x(:,tp(1):tp(2));
y = y(:,tp(1):tp(2));

%% raw
col = fp.c_area(idx==0,:);
PlotActivity(x,area_val(idx==0),area_label(idx==0),col)

%% reorganize by the first dimension
BB = B(:,2);
if sum(BB>0)<sum(BB<0)
    BB = -1*BB; 
end
[~,nidx] = sort(BB,'descend');
PlotActivity(x(nidx,:),area_val(idx==0),area_label(idx==0),col)
% add a line over the top showing average activity

%% 
figure; hold on; 
y = x'*B(:,1:3);
y = movmean(y,15,1);
plot(y(:,1),y(:,2));
figure; hold on; 
for i = 1:size(y,1)
    plot(y(i,1),y(i,2),'.');
    pause(0.1);
end


end %function end



%% subroutines
function PlotActivity(x,area_val,area_label,col)
fp = fig_params_cortdynamics;
figure; hold on;
%%uncomment for imagesc method
imagesc(x,[0 1]); colorbar; 
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
    plot([-5 -5],idx,'color',col(i,:),'linewidth',3); 
    text(-10, idx(1)+(idx(2)-idx(1))/2, area_label{i},'FontWeight','bold','HorizontalAlignment','right','Color',fp.c_area(i,:),'fontsize',fp.font_size,'FontName',fp.font_name)
    
end
set(gca,'YAxisLocation','right');
set(gca,'xtick',0:(15*10):size(x,2),'xticklabel',(0:(15*10):size(x,2))/15)
xlabel('time (seconds)');
ylabel('neurons')
fp.FormatAxes(gca); box on
set(gca,'Clipping','off')
end