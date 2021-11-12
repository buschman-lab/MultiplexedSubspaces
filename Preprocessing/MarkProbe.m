function MarkProbe(opts_fn,probe_label)
%Camden MacDowell - timeless
%Manual marking of probe surface trajectory
%opts_fn is a cell array 

if nargin <2
    probe_label = {'1 Par','2 SS','3 Vis','4 PFC'};
end

for cur_fn = 1:numel(opts_fn) %file loop
    opts = load(opts_fn{cur_fn},'opts','prepro_log');
    if isfield(opts,'opts') %will either be prepro or opts depending on where in pipeline
        opts = opts.opts;  
    elseif isfield(opts,'prepro_log')
        opts = opts.prepro_log; 
    else
        error('no options or prepro_log'); 
    end
    
    while 1 %redo loop
        dlg = 'Yes';       
        position = [];
        figure('name','Select two points along probe IN ORDER. First at insertion','units','normalized','position',[0 0 1 1]);
        imagesc(opts.cropped_alligned_img); colormap gray; hold on;
        probe_coords = {}; 
        count=1;
        while 1 %probe loop            
            switch dlg
                case 'Yes' %mask additional region
                    if ~isempty(position) %mark the previous probes
                        plot(position(:,1),position(:,2),'color','g','linewidth',2)
                        plot(position(1,1),position(1,2),'marker','.','markersize',20,'color','r','linestyle',':')      
                        text(position(1,1),position(1,2),probe_label{count},'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',10,'color','b');
                        count = count+1;
                    end                                                              
                    midline = imline; %drawline functionality in 2018b isn't quite there yet so don't use
                    position = wait(midline); 
                    probe_coords = cat(1,probe_coords,{round(position/opts.spatial_bin_factor)});                    
               case 'No'; break %end dlg
            end
            choice = questdlg(sprintf('Mark additional probes (in experiment order)'),...
                'mark ADDITIONAL probes',...
                'Yes','No','No');
            switch choice; case 'Yes';case 'No'; break; end        
        end  %probe while       
        plot(position(:,1),position(:,2),'color','g','linewidth',2)
        plot(position(1,1),position(1,2),'marker','.','markersize',20,'color','r','linestyle',':')
        text(position(1,1),position(1,2),probe_label{count},'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',10,'color','b');
        choice = questdlg(sprintf('Redo?'),...
            'mark ADDITIONAL probes',...
            'Yes','No','No');
        switch choice; case 'Yes';case 'No'; break; end    
    end %redo loop
    close;     
    
    %append the new info
    if ~isempty(probe_coords)
        save(opts_fn{cur_fn},'probe_coords','-append');
    end
    
end %file loop

end %function

    