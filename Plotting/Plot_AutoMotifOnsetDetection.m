function Plot_AutoMotifOnsetDetection(savedir)
%Camden MacDowell - timeless
%see MotifOnset for details. 
%takes ~1 hour to run

if nargin <1; savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\MotifOnset'; end
if~exist(savedir,'dir'); mkdir(savedir); end

fp = fig_params;

%Match motif fits to the grab data
[motif_fits,~] = GrabFiles('\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
[orig_data,~] = GrabFiles('\w*processed\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging'});


%preallocate
onsets = cell(numel(orig_data),44); %44 is the number of chunks in 2021 data
onsets_mask = cell(numel(orig_data),44); %44 is the number of chunks in 2021 data
thresh_level = cell(numel(orig_data),44);
h_wf = cell(numel(orig_data),44);
rho = cell(numel(orig_data),44);
rho_null = cell(numel(orig_data),44);

%loop through original data
for cur_file = 1:numel(orig_data)
    data = load(orig_data{cur_file},'data_test','data_train');    
    [~,fn] = fileparts(orig_data{cur_file});
    %loop through each fit
    idx = find(contains(motif_fits,fn)==1);
    for cur_fit = 1:numel(idx)
       fit_data = load(motif_fits{idx(cur_fit)},'w','H');
        
       %match to test or train data
       [~,temp_fn]= fileparts(motif_fits{idx(cur_fit)});
       
       if contains(temp_fn,'test')
           chunk = extractBetween(temp_fn,'chunk_','test');
           [onsets{cur_file,cur_fit},onsets_mask{cur_file,cur_fit},thresh_level{cur_file,cur_fit},h_wf{cur_file,cur_fit},rho{cur_file,cur_fit},rho_null{cur_file,cur_fit}]=MotifOnset(fit_data.w,fit_data.H,data.data_test(:,:,str2double(chunk{1})));
       elseif contains(temp_fn,'train')
           chunk = extractBetween(temp_fn,'chunk_','train');           
           [onsets{cur_file,cur_fit},onsets_mask{cur_file,cur_fit},thresh_level{cur_file,cur_fit},h_wf{cur_file,cur_fit},rho{cur_file,cur_fit},rho_null{cur_file,cur_fit}]=MotifOnset(fit_data.w,fit_data.H,data.data_train(:,:,str2double(chunk{1})));
       else
           error('no test or train data in file');
       end
    end %fit loop
    
end %file loop

%concatenate and plot by motif
n = size(rho{1},1);
motif_idx = repmat(1:n,1,numel(rho));

%rho
rho_cat = cat(1,rho{:});

%rho_null
rho_null_cat = cat(1,rho_null{:});

%thresh level
thresh_level_cat = cat(2,thresh_level{:});

%h waveform
h_wf_cat = cellfun(@(x) cell2mat(x)',h_wf,'UniformOutput',0);
h_wf_cat = cat(1,h_wf_cat{:});


%figures
for i = 1:n
    close all
    figure; hold on; 
    line([1 size(rho_cat,2)],[0 0],'linewidth',1,'color','k');
    plot(rho_null_cat(motif_idx==i,:)','color',[0.5 0.5 0.5 0.15],'linewidth',0.5)    
    plot(rho_cat(motif_idx==i,:)','color',[0.8500 0.3723 0.0078 0.15],'linewidth',0.5)    
    xlabel('Threshold Level');
    ylabel('rho');
    set(gca,'xlim',[1,size(rho_cat,2)])
    fp.SetTitle(gca,sprintf('Adaptive threshold Motif %d',i));
    fp.FormatAxes(gca)
    fp.FigureSizing(gcf,[3 3 6 4],[])
    
    %plot the thresholds used.
    figure; hold on;
    histogram(thresh_level_cat(motif_idx==i),'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none');
    xlabel('Threshold Level');
    ylabel('# Epochs');
    set(gca,'xlim',[1,size(rho_cat,2)])
    fp.SetTitle(gca,sprintf('Threshold Level Motif %d',i));
    fp.FormatAxes(gca)    
    fp.FigureSizing(gcf,[3 3 3 4],[])

    %plot the H timecourse
	figure; hold on;     
    plot(h_wf_cat(motif_idx==i,:)','color',[0.5 0.5 0.5 0.15],'linewidth',0.5)    
    xlabel('Time (frames)');
    ylabel('H Weight');
    set(gca,'xlim',[1,size(h_wf_cat,2)])
    fp.SetTitle(gca,sprintf('H Waveform Motif %d',i));
    fp.FormatAxes(gca)
    fp.FigureSizing(gcf,[3 3 6 4],[])
        
    %save off
    saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('AdaptiveThresholding_motif%d',i),savedir,0); close all
end    










