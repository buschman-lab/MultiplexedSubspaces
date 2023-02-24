function PlotExampleImagingTrace(cur_rec,tp)
fp = fig_params_cortdynamics;
%camden macdowell
if nargin <2; tp = [100,100+(60*15*1)]; end

%VANESSA added to fix pathing please feel free to delete later
if all(getenv('username') == 'roser')
    [rec_name,ImgPath,ImgProbeLoc,~,motif_fits] = LoadDataDirectories_VR(cur_rec);
else
    [rec_name,ImgPath,ImgProbeLoc,~,motif_fits] = LoadDataDirectories(cur_rec);
end
%load the deconvolve dff form the probe insertion sites
[dff_probe,offset,dff_contra] = LoadInsertionSiteDFF(ImgPath,ImgProbeLoc);%probe offset is a few pixel offset in the probe insertion location used for deconvolution to adjust for the fact that the probes are inserted at an angle so the first recorded cell bodies are m/l and a/p shifted from the exact probe insertion site at the surface. 

%load the probe coordinates 
probe_loc = load(ImgProbeLoc); probe_loc = probe_loc.probe_coords;

x = dff_probe(tp(1):tp(2),:);
x = x-min(x);
y = dff_contra(tp(1):tp(2),:);
y = y-min(y);
figure; hold on; 
%hardcode offsets
hemi_offset = [3,11,19,27];
loc_offset = [0,8,16,24];
for i = 1:size(dff_probe,2)
    temp = 1:size(x,1);
    plot(temp,x(:,i)+loc_offset(i),'linestyle','-','linewidth',2,'color',[0.75 0.75 0.75])
    plot(temp,y(:,i)+hemi_offset(i),'linestyle','-','linewidth',2,'color',[0.75 0.75 0.75])
end

xlim([-5,size(x,1)+0.5])
set(gca,'xtick',0:(15*10):(tp(2)-tp(1)),'xticklabel',(0:(15*10):(tp(2)-tp(1)))/15)
xlabel('time (seconds)');
ylabel({'Estimated FR','(deconvolved \DeltaF/F\sigma)'})
fp.FormatAxes(gca); box off



end %function end
