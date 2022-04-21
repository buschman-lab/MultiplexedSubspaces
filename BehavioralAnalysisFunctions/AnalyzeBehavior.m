function [me,roi_names,w,U,xx,yy] = AnalyzeBehavior(cur_rec,svdFlag)
if nargin <2; svdFlag =0; end
%Camden MacDowell - timeless

%Load all chunks for a recording
rec_name = LoadDataDirectories(cur_rec);
if ispc
    [fn,~] = GrabFiles([rec_name,'\w*pro_chunk\w*.mat'],0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BehavioralData\'});
else
    [fn,~] = GrabFiles([rec_name,'\w*pro_chunk\w*.mat'],0,{'/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/BehavioralData/'});
end

if svdFlag==1
    data = cellfun(@(x) load(x,'roi_names','xx','yy','u','me','face'),fn,'UniformOutput',0);

    %pan-recording SVD
    u = cellfun(@(x) x.u,data,'UniformOutput',0);
    u = cat(2,u{:});

    U = svdecon(u);
    %get the explained variance of the dimensions
    U = normc(u); %normalize by column

    %compile face motion video (to get temporal weightings of eigen faces)
    vid = cellfun(@(x) x.face,data,'UniformOutput',0);
    vid = cat(2,vid{:});

    % get the temporal weightings
    w = U'*vid; 
    w(1:1000,:);
    U = U(:,1:1000);    
elseif svdFlag==2 %perform PCA (more interpretable IMO)
    data = cellfun(@(x) load(x,'roi_names','xx','yy','u','me','face'),fn,'UniformOutput',0);

    %compile face motion video (to get temporal weightings of eigen faces)
    vid = cellfun(@(x) x.face,data,'UniformOutput',0);
    vid = cat(2,vid{:});
    w = vid; 
    U = [];
else
   data = cellfun(@(x) load(x,'roi_names','me','xx','yy'),fn,'UniformOutput',0);
    w = [];
    U = [];        
end

%compile the motion energy 
me = cellfun(@(x) x.me,data,'UniformOutput',0);
me = cat(2,me{:});
roi_names = data{1}.roi_names;     

xx = data{1}.xx;
yy = data{1}.yy;

end

%uncomment to watch it (gut check)
% xx = data{1}.xx;
% yy = data{1}.yy;
% vid = reshape(vid,xx,yy,size(vid,2));
% for i = 1:size(vid,3)
%    imagesc(vid(:,:,i)); colormap gray; 
%    title(sprintf('%d',i));
%    pause(0.1);
% end