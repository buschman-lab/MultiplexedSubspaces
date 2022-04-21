function PlotBasisMotifs(savedir,cur_mouse)
%Camden 2022
if nargin <2; cur_mouse=[]; end

if isempty(cur_mouse)
    W_basis = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\MotifDiscovery\Mouse_basis_motifs.mat','W_basis');   
else
    W_basis = load(['Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\MotifDiscovery\','Mousepermouse_basis_motifs',num2str(cur_mouse),'.mat'],'W_basis');   
end
W_basis = W_basis.W_basis;
fp = fig_params_cortdynamics;
W_basis(:,fp.noisemotif,:) = [];

%scale for each, since each has a diff number of pixels
cmax = 98.5*ones(size(W_basis,2),1);
cmax(14) = 97.5;
cmax(8) = 97.5;
cmax(6) = 98.;

for i = 1:size(W_basis,2)
    temp = squeeze(W_basis(:,i,:));
    temp = reshape(temp,[68 68 size(temp,2)]);
    temp(temp==0)=NaN;
    figure;
    cvals = [0,prctile(temp(:),cmax(i))];
    for j = 1:size(temp,3)
        PlotMesoFrame(temp(:,:,j),'caxis',cvals,'caxis_flag',1);
        title(sprintf('motif %d | %d',i,j));
%         pause()
        if isempty(cur_mouse)
            saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_frame%d',i,j),savedir,0); close all
%             saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_frame_sided%d',i,j),savedir,0); close all
        else
            saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('%dmotif%d_frame%d',cur_mouse,i,j),savedir,0); close all
        end
    end
    
end
    


end