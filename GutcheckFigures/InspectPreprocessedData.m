function InspectPreprocessedData(filepath,type)
%camden macdowell 2020

%%
switch type
    case 'preprocessed'
        load(filepath,'dff');
        figure;
        for j = 10000:size(dff,3) %I like to start in the middle to make sure nothing went funky during the allignment and hemodynamic correction
          imagesc(dff(:,:,j),[0 3]);      
          axis equal; axis off
          colormap magma
          title(sprintf('frame %d',j),'fontsize',12,'fontweight','normal');
          pause(0.02);
        end
        
    case 'postprocessed'
        load(filepath,'data_norm','nanpxs');
        dff = conditionDffMat(data_norm',nanpxs);
        figure;
        for j = 20000:size(dff,3) %I like to start in the middle to make sure nothing went funky during the allignment and hemodynamic correction
          imagesc(dff(:,:,j),[0 .1]);      
          axis equal; axis off
          colormap magma
          title(sprintf('frame %d',j),'fontsize',12,'fontweight','normal');
          pause(0.02);
        end
        
end
        
%%