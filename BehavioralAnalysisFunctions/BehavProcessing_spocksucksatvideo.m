function BehavProcessing_spocksucksatvideo(cur_chunk,cur_rec,fn)
%Camden MacDowell - timeless
% Follows Stringer et a., 2019 | Science

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
end

[rec_name,~,~,~,~,FaceCam,~,~] = LoadDataDirectories(cur_rec);

ndim = 100;

fprintf('Working on chunk %d',cur_chunk);
%load the parse video data
load(fn,'chunk','roi','roi_names');
t = chunk{cur_chunk};

v = VideoReader(FaceCam);
pause(10);
v.NumFrames
fprintf('current video time %d',v.CurrentTime);

% load the chunk since SPOCK can't use read (takes longer and longer
% unforunately for later sections
CurFrame = 0;
cur_t=1;
data = NaN(v.Height*v.Width,numel(t));
while 1 %hasFrame(v)    
    fprintf('\n cur frame %d',CurFrame);
    CurImage = readFrame(v);
    CurFrame = CurFrame+1;
    if ismember(CurFrame, t)
        fprintf('\n aquired frame %d',CurFrame);
        temp = rgb2gray(CurImage);
        data(:,cur_t) = temp(:);
        if cur_t==size(data,2)
            break
        end
        cur_t = cur_t+1;        
    end
end
fprintf('Done loading video');
xx = roi{1}.position(4)+1;
yy = roi{1}.position(3)+1;
zz = numel(t)-1;

fprintf('\nParsing Video...')
me = NaN(numel(roi)-1,numel(t)-1);
face = NaN(xx*yy,numel(t)-1);
for i = 1:size(data,2)-1
%     img = abs(double(rgb2gray(read(v,t(cur_t))))-double(rgb2gray(read(v,t(cur_t+1)))));
    img = reshape(data(:,i),v.Height,v.Width)-reshape(data(:,i+1),v.Height,v.Width);
    %Compute summed ME individual rois 
    for cur_roi = 2:numel(roi)
        me(cur_roi,i)=nansum(imcrop(img,roi{cur_roi}.position),'all');
    end     
    
    %compute pca of the face
    a = imcrop(img,roi{1}.position);    
    face(:,i) = a(:);
end
fprintf('\nDone Parsing Video \n saving now...')

%SVD hangs on spock... not sure why, it should take like 1 second. So save
%off face instead and just stop here. 
% fprintf('\nDone Parsing \nVideo Computing SVD with econ...')
% %do svc and get the usv
% tic; [u,~,~] = svdecon(face); toc
% u = normc(u); %normalize by column
% fprintf('\nSaved - Done with SVD');

savedir = fileparts(fn);
f = [savedir,filesep,rec_name,sprintf('pro_chunk%d.mat',cur_chunk)];
save(f,'face','me','cur_rec','ndim','zz','xx','yy','roi_names','roi','CurFrame','t');

fprintf('\nSaved - Done');

end


