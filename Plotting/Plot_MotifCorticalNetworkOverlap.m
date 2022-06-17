function Plot_MotifCorticalNetworkOverlap(data)

%for a brain region and recording, load the maps for a motif
cur_rec = 1;
cur_a =1; 

%%
overlap = cell(6,8);
idx = cell(6,8);
d = cell(6,8);
for cur_rec = 1:6
    fprintf('\nWorking on %d',cur_rec);
    for cur_a = 1:8
        map = cat(2,data{cur_rec});
        index = [map.cur_a]==cur_a;
        map = map(index);
        y = cell(1,numel(map));
        for i = 1:numel(map)
            y{i} = getSig(map(i));
        end
        [overlap{cur_rec,cur_a},a,b] = CompareMotifs(y);
        a(eye(size(a))==1)=NaN;
        b(eye(size(b))==1)=NaN;
        idx{cur_rec,cur_a} = a;
        d{cur_rec,cur_a} = b;
    end
end
x = cat(1,idx{:});
histogram(x(:))
nanmean(x(:))
bootci(1000,@nanmean,x(:))

y = cat(1,d{:});
bootci(1000,@nanmean,y(:))
histogram(y(:))
nanmean(y(:))

x = cat(1,overlap{~cellfun(@isempty,overlap)});
x(isnan(x))=[];
close all;histogram(x(:))
% bootci(1000,@nanmean,x(:))

%%
figure; hold on; 
x = cat(1,overlap{~cellfun(@isempty,overlap)});
x(isnan(x))=[];
histogram(x(:)*100,'facecolor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5,'binwidth',100*0.05)
xlim([0 100])
xlabel('% Overlap')
ylabel({'cross-motif','comparisons'})
fp.FormatAxes(gca); box on; 
fp.FigureSizing(gcf,[3 2 2.5 2.5],[2 10 10 10])
[~,stats] = pairedBootstrap(x(:)*100,@nanmean);
title(sprintf('%0.2f %0.2f %0.2f',stats.mean,stats.ci(1),stats.ci(2)),'Fontweight','normal');

%% plot across dimensions
% x = cat(4,overlap{~cellfun(@isempty,overlap)});
% x = permute(x,[1,2,4,3]);
% x = reshape(x,size(x,1)*size(x,2)*size(x,3),size(x,4));
% xboot = pairedBootstrap(x,@nanmean);
% figure; hold on; 
% CompareViolins(xboot',fp,'plotspread',0,'divfactor',0.2,'plotaverage',1,'distWidth',0.75);

%plot the rearangement

end


function y = getSig(x)
    y = NaN(68,68,10);
    for i = 1:10
       temp = x.rho_all(:,:,i);
       temp(abs(temp)<x.sig_thresh(i))=0;
       temp(abs(temp)>0)=1;  
       if sum(temp(:)==1)<(0.01*sum(~isnan(temp(:)))) %ignore if less than 1% of pixels (~<20)
           temp = NaN;
       end
       y(:,:,i) = temp;      
    end
    temp = reshape(y,size(y,1)*size(y,2),size(y,3));
    bad_idx = sum(isnan(temp))==size(temp,1);
    y(:,:,bad_idx)=[];
    
end


function [overlap,idx,d] = CompareMotifs(y)
    overlap = NaN(numel(y),numel(y),10);
    idx = NaN(numel(y),numel(y));
    d = NaN(numel(y),numel(y));
    for i = 1:numel(y)
        for j = 1:numel(y)   
            if i>j
                [a,idx(i,j),d(i,j)] = ComputeOverLap(y{i},y{j});
%                 a(eye(size(a))==1) = nan;
                overlap(i,j,:) = a;
            end
        end
    end
    
end

function [overlap,idx,d] = ComputeOverLap(x,y)
   x = reshape(x,size(x,1)*size(x,2),size(x,3));
   y = reshape(y,size(y,1)*size(y,2),size(y,3));
   overlap = NaN(1,10);
   idx = NaN(1,10);
   for i = 1:size(x,2)
      temp = arrayfun(@(n) (2*sum(x(:,i)+y(:,n)==2))/(sum(x(:,i)==1) + sum(y(:,n)==1)),1:size(y,2),'UniformOutput',1);
      [overlap(i),idx(i)] = max(temp);
   end
   overlap = nanmean(overlap);
   %get the distance
   d = nanmean(abs(idx(1:end)-(1:10)));
   idx = corr(idx(2:end)',[2:10]','type','Kendall');
end




% 
% 
% map = cat(2,data{:});
% %go through each one
% x = cell(1,size(map,2));
% rho= cell(1,size(map,2));
% ov = cell(1,size(map,2));
% for j = 1:size(map,2)
%     y = map(j).rho_all;
%     p = map(j).sig_thresh;    
%     ysig = y;
%     for i = 1:10
%        temp = ysig(:,:,i);
%        temp(abs(temp)<p(i))=0;
%        temp(abs(temp)>0)=1;
%        ysig(:,:,i) = temp;
%     end
%     
%     %flip rho_all so strongest is positive
%     for i = 1:size(y,3)
%         if nansum(y(:,:,i),'all') < nansum(-1*y(:,:,i),'all')
%             y(:,:,i)=-1*y(:,:,i);
%         end
%     end    
%     % get the overlap between all combineations
%     ysig = reshape(ysig,size(ysig,1)*size(ysig,2),size(ysig,3));
%     ytemp = reshape(y,size(y,1)*size(y,2),size(y,3));
%     x{j} = NaN(10,10);
%     rho{j} = NaN(10,10);
%     temp = NaN(size(ysig,1),100);
%     %get the overlap
%     COUNT=1;
%     for i = 1:size(ysig,2)
%         for ii = 1:size(ysig,2)
%             if ii>i
%                 x{j}(i,ii) = (2*sum((ysig(:,i)+ysig(:,ii))==2))/(sum(ysig(:,ii)==1)+sum(ysig(:,i)==1));
%                 rho{j}(i,ii) = corr(ytemp(:,i),ytemp(:,ii),'rows','complete');
%                 temp(:,COUNT) = (ysig(:,i)+ysig(:,ii))==2;
%                 COUNT = COUNT+1;
%             end
%         end
%     end  
%     ov{j} = nansum(temp,2);
% end