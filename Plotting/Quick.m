%some basic files checks on imaging

fn = GrabFiles('dff_combined.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging'});

%confirm nothing funky with frames
t = NaN(1,numel(fn));
for i = 1:numel(fn)
    matObj = matfile(fn{i});
    t(i) = size(matObj,'dff',3);
end

N = t(1);
%plot some example data
figure;
for cur_fn = 1:numel(fn)  
    matObj = matfile(fn{cur_fn});
    for i = 2500:3000
        imagesc(matObj.dff(:,:,i),[-4 4]);
        colormap magma
        title(sprintf('%d %d',cur_fn,i));
        pause(0.01);
    end
end
    

