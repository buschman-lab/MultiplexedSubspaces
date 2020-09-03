function MotifToGif(W,save_path,varargin)
%Camden - timeless
%W is a px x px x L so this also works from raw data movies

opts.limit = 95; %percentile of limit
opts.delay = 0.1; %delay between frames
opts.type = 'simple'; %how fancy should we be?
opts = ParseOptionalInputs(opts,varargin); 

switch opts.type
    case 'simple'
        figure; 
        im = cell(1,size(W,3));
        for i = 1:size(W,3)
            cla
            imagesc(W(:,:,i),[0 prctile(W(:),opts.limit)]);
            colormap magma; axis equal;
            im{i} = frame2im(getframe(gcf));            
        end
    case 'fancy'
        %see MakeRepitoireFig in the original code repo for first paper
    otherwise 
        error('unknown type');
end

for idx = 1:numel(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,save_path,'gif','LoopCount',Inf,'DelayTime',.5);
    else
        imwrite(A,map,save_path,'gif','WriteMode','append','DelayTime',opts.delay);
    end
end
close(gcf)
end %end function