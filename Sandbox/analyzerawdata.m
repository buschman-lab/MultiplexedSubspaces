x = (1:size(features_downsampled,1));
figure('position',[ 358         152        1162         826]); hold on
b = subplot(211); hold on; 
plot(x/(60*13),features_downsampled);
title('behavioral states')
setFigureDefaults
set(gca,'xgrid','on','GridColor','k')
legend('location','eastoutside')
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) 20 6])

a = subplot(212); hold on; 
plot(x/(60*13),H_weight);
xlabel('Time (min)')
title('motifs')
setFigureDefaults
set(gca,'xgrid','on','GridColor','k')
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) 20 6])

linkaxes([a,b],'x');


%%
n_shuf=1000;
rng('default')
data = load([fn_path fn_widefield],'w','data_test','H');
w = data.w(:,expvaridx,:); 
H = data.H(expvaridx,:);

H_weight = NaN(size(H))';
for cur_motif = 1:size(H,1)
    H_weight(:,cur_motif) = helper.reconstruct(nanmean(w(:,cur_motif,:),1),H(cur_motif,:)); 
end %motif loop

%%
H_weight = H'; 

data_test = load(fn_widefield,'data_test');

H_weight = nanmean(data_test.data_test,1)';
%%
temp = [];
for i = 1:size(H_weight,2)
    temp(:,i) = H_weight(:,i)-nanmean(H_weight(:,i));
end
temp_behav = features_downsampled-nanmean(features_downsampled);

num_features=3;
% Plot the autocorrelation of factors 
figure('units','normalized','position',[0 0 1 1]); hold on;
for cur_motif = 1:size(H_weight,2)
    subplot(2,7,cur_motif); hold on
    for i = 1:num_features
        %plot pdf 
        [xc,lags] = xcorr(temp_behav(:,i),temp(:,cur_motif),1000*13,'coeff'); 
%         idx = [ceil(numel(lags)/2):numel(lags)];
        plot(lags/13,xc,'linewidth',2,'color',col(i,:));   
        title(sprintf('Motif %d',cur_motif),'FontName','Arial','FontSize',16,'FontWeight','normal')    
        ylim([min(xc) 0.25]);
%         xlim([0 max(lags/13)])
        ylabel('Rho')
        xlabel('time (s)')

       setFigureDefaults;       
    end
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) 3 3])
end












