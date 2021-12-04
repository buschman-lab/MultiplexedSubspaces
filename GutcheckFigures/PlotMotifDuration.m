function PlotMotifDuration(file_list,parameter_class)
%Camden MacDowell - timeless
%loads the cost and regularization and averages across chunks

if nargin <1
   [file_list,~] = GrabFiles('\w*fit\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\MotifDiscovery'}); %select the preprocessed data (not the '_processed');
end

if nargin <2; parameter_class = 'general_params_corticaldynamics'; end 
    
gp = loadobj(feval(parameter_class)); 

data = cellfun(@(x) load(x,'w'),file_list,'UniformOutput',0);

[~,~,w] = cellfun(@(x) ShiftW(x.w), data,'UniformOutput',0);
w = cellfun(@(x) squeeze(nanmean(x,1)), w,'UniformOutput',0);
w = cat(1,w{:});

%pad and shift COM to middle
fp = fig_params;
figure; hold on; 
shadedErrorBar(1:size(w,2),nanmean(w),sem(w),'lineprops',{'color','r'})
set(gca,'xscale','linear');
xlabel('frames'); ylabel('activity');
title('Motif duration','Fontweight',fp.font_weight,'Fontsize',fp.font_size)


end