function ImpactOfNormalization()
fp = fig_params_cortdynamics;

nval = 1000;
fr = linspace(0,5,nval);
a = NaN(nval,nval);
b = NaN(nval,nval);
for i = 1:nval
    for j = 1:nval
       a(i,j) = fr(i)/(fr(j)+1); 
       b(i,j) = fr(i)-fr(j); 
    end
end

%divisive normalization
figure; hold on; 
imagesc(a,[0 3]); c = colorbar; colormap magma
ylabel(c,'Normalized FR');
set(gca,'XScale','log','YScale','log')
tickval = get(gca,'xtick');
set(gca,'xtick',tickval,'xticklabel',round(fr(tickval),2),'XTickLabelRotation',45);
set(gca,'ytick',tickval,'yticklabel',round(fr(tickval),2));
ylabel({'Trial FR','(spikes/sec)'}); 
xlabel({'Baseline FR','(spikes/sec)'}); 
title({'Divisive normalization increases','impact of sparse an selective neurons'},'fontweight','normal');
fp.FormatAxes(gca);  box on; 
fp.FigureSizing(gcf,[3 3 2 2],[10 10 10 12])


%subtraction normalization
figure; hold on; 
imagesc(b,[-5 5]); c = colorbar; colormap magma
ylabel(c,'Normalized FR');
set(gca,'XScale','log','YScale','log')
tickval = get(gca,'xtick');
set(gca,'xtick',tickval,'xticklabel',round(fr(tickval),2),'XTickLabelRotation',45);
set(gca,'ytick',tickval,'yticklabel',round(fr(tickval),2));
ylabel({'Trial FR','(spikes/sec)'}); 
xlabel({'Baseline FR','(spikes/sec)'}); 
title('Subtraction normalization','fontweight','normal');
fp.FormatAxes(gca);  box on; 
fp.FigureSizing(gcf,[3 3 2 2],[10 10 10 12])

end