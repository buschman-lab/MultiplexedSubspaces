function Figure1i(Dshared,Dlocal)
fp = fig_params_cortdynamics;
col = arrayfun(@(n) fp.c_area(n,:),1:size(fp.c_area,1),'UniformOutput',0);

LvsIR(Dshared,Dlocal,col); 
fp.FormatAxes(gca);
box on; grid off
fp.FigureSizing(gcf,[3 2 3.75 3.75],[10 10 15 10]) 
xlim([1 7]); ylim([1 7])

end

function fh = LvsIR(X,Y,col)
fh = figure; hold on; 
for i = [1,2,4,5,6,7,8,3] %reorder so that prelimbic is plotted at end since yellow and gets covered
    x = squeeze(X(:,:,i));
    x = reshape(x,size(x,1)*size(x,2),size(x,3));
    y = squeeze(Y(:,:,i));
    y = reshape(y,size(y,1)*size(y,2),size(y,3));
    
    %bootstrap the distribution
    [x,~] = pairedBootstrap(x,@nanmean); 
    [y,~] = pairedBootstrap(y,@nanmean); 
    
    %plot the cloud
    s1 = scatter(x(:),y(:),'o','MarkerEdgeAlpha',0.1);  
    s1.SizeData = s1.SizeData/10; 
    s1.MarkerEdgeColor = col{i};   

    ylabel('Local Dimensions');
    xlabel({'Shared Subspace','Dimensions'});    
    plot([0 15],[0 15],'linestyle',':','color',[0.8 0 0],'linewidth',2)
    plot([0:15],2*[0:15],'linestyle',':','color',[0.5 0.5 0.5],'linewidth',2)      
end
    
end
