function BehavProcessing(cur_rec,fn)
%Camden MacDowell - timeless
% Follows Stringer et a., 2019 | Science
% Videos are big so loads and processes in chunks

[rec_name,~,~,~,~,FaceCam,~,~] = LoadDataDirectories(cur_rec);

v = VideoReader(FaceCam);
ndim = 100;

load(fn,'chunk','roi','roi_names');
CurFrame = 1;
settime = 0; 
for cur_chunk = 1:numel(chunk)
    tic
    fprintf('Working on chunk %d',cur_chunk);
    %load the parse video data
    t = chunk{cur_chunk};

    %load the video, starting where you had previously stopped    
    [CurFrame, data] = LoadVideo(v,CurFrame-1,t,settime);
    settime = v.CurrentTime; 
    
    xx = roi{1}.position(4)+1;
    yy = roi{1}.position(3)+1;
    zz = numel(t)-1;

    fprintf('\nParsing Video...')
    me = NaN(numel(roi)-1,numel(t)-1);
    face = NaN(xx*yy,numel(t)-1);
    for i = 1:size(data,2)-1
        img = reshape(data(:,i),v.Height,v.Width)-reshape(data(:,i+1),v.Height,v.Width);
        %Compute summed ME individual rois 
        for cur_roi = 2:numel(roi)
            me(cur_roi,i)=nansum(imcrop(img,roi{cur_roi}.position),'all');
        end     

        %compute pca of the face
        a = imcrop(img,roi{1}.position);    
        face(:,i) = a(:);
    end

    %SVD hangs on spock... not sure why, it should take like 1 second. So save
    %off face instead and just stop here. 
    fprintf('\nDone Parsing \nVideo Computing SVD with econ...')
    %do svc and get the usv
    [u,~,~] = svdecon(face);
    u = normc(u(:,1:ndim)); %normalize by column
    fprintf('\nDone with SVD');

    savedir = fileparts(fn);
    f = [savedir,filesep,rec_name,sprintf('pro_chunk%d.mat',cur_chunk)];
    save(f,'u','face','me','cur_rec','ndim','zz','xx','yy','roi_names','roi','CurFrame','t');

    fprintf('\nSaved - Done');
    tic
end

end

function [CurFrame, data] = LoadVideo(v,CurFrame,t,settime)
cur_t=1;
v.CurrentTime = settime; 
data = NaN(v.Height*v.Width,numel(t));
while hasFrame(v)    
    CurImage = readFrame(v);
    CurFrame = CurFrame+1;
    if ismember(CurFrame, t)
        temp = rgb2gray(CurImage);
        data(:,cur_t) = temp(:);
        if cur_t==size(data,2)
            break
        end
        cur_t = cur_t+1;        
    end
end
end
