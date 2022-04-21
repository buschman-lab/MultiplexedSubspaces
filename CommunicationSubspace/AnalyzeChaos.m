function AnalyzeChaos(folder,rec_name)

folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\chaos';
rec_name = 'Mouse 331 Recording 1';
[fn,~] = GrabFiles([rec_name '\w*Chaos\w*.mat'],0,{folder}); 
c = cellfun(@(x) load(x,'c'),fn);
label = load(fn{1},'area_label'); label = label.area_label;
c_all = cat(1,c(:).c);

rat_chaos = NaN(size(c_all));
for i = 1:size(c_all,1)
    for j = 1:size(c_all,2)
       temp = cellfun(@(x) x.result,c_all{i,j},'UniformOutput',false);
       rat_chaos(i,j) = sum(ismember(temp,'chaotic'))/numel(temp);
    end
end

figure; hold on; 
plot(rat_chaos*100)
legend(label);
ylabel('% chaotic neurons');
xlabel('time (min)');
xval = round(linspace(0,(40000*2/15/60),39));
set(gca,'XTick',1:5:40,'XTickLabel',xval(1:5:40));
title('absolute chaos','fontweight','normal')

%correlation in chaos between areas
figure; hold on; 
rho = corr(rat_chaos-nanmean(rat_chaos));
imagesc(rho,[-0.1 0.5]); colorbar; colormap magma
set(gca,'XTick',[1:8],'XTickLabel',label)
set(gca,'YTick',[1:8],'YTickLabel',label)
title('correlation in chaos')

%average chaos versus average dimensionality
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';
[fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_mo\w*.mat'],0,{folder}); 
data_rrr = cellfun(@(x) load(x),fn);
[area_label,d_opt,d,n] = localDimensionality(data_rrr); 

%correlation between chaos and dimensionality
figure; hold on; 
d = nanmedian(d_opt,2);
lm = fitlm(nanmedian(rat_chaos),d');
plot(lm)
rho = corr(d,nanmedian(rat_chaos)');
xlabel('average Chaos %'); ylabel('dimensionality')
 title({'relationship between chaos and dimensionality',sprintf('R %0.2f',rho)})

%correlation in chaos and number of neurons
figure; hold on; 
n = nanmedian(n,2);
lm = fitlm(nanmedian(rat_chaos),n');
plot(lm)
rho = corr(n,nanmedian(rat_chaos)');
xlabel('average Chaos %'); ylabel('# nuerons')
title({'relationship between chaos and # neurons',sprintf('R %0.2f',rho)})




end %funciton end












