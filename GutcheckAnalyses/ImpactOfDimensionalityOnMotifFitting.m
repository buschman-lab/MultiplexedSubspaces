function ImpactOfDimensionalityOnMotifFitting(data_train,data_test)
%Camden MacDowell - timeless
%Looks at the impact the dimensions of input matrix X on number of
%discovered motifs and computational time

%test impact of duration
dur = round(15*60*[0.25, 0.5, 1, 2, 3, 4, 5]);
compt_time = NaN(1,numel(dur));
dur_pev = NaN(numel(dur),2);
for i = 1:numel(dur)
    fprintf('\n\t Working on Duration %d of %d',i,numel(dur));
    data_train_temp = data_train(:,1:dur(i));
    data_test_temp = data_test(:,1:dur(i));
    tic
    [stats_test, stats_train, w, h, gp, fh] = FitMotifs(data_train_temp,data_test_temp,'maxiter',50,'lambda_range', sort(logspace(-1,-4,5), 'ascend'));
    t = toc;
    dur_pev(i,1) = stats_train.pev;
    dur_pev(i,2) = stats_test.pev;
    compt_time(i) = t; 
end

figure; hold on; 
plot((dur_pev(:,1)-min(dur_pev(:,1)))/(max(dur_pev(:,1))-min(dur_pev(:,1))),'r')
plot((dur_pev(:,2)-min(dur_pev(:,2)))/(max(dur_pev(:,2))-min(dur_pev(:,2))),'b')
plot((compt_time-min(compt_time))/(max(compt_time)-min(compt_time)),'k')





%test impact of number of pixels
