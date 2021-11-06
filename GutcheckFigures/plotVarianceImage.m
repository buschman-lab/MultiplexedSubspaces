function plotVarianceImage(stack,var_type,probe_locations)
%Camden MacDowell - timeless
%Makes a nice plot of the variance over time across images
%stack is the processed recording state
%if given probe locations, will plot those
if nargin <2; var_type='std'; end
if nargin <3; probe_locations = []; end

stack(isnan(stack))=0;
switch var_type
    case 'std'
        frame = imgaussfilt(nanvar(stack,[],3),'filtersize',[5,5]);        
    case 'var'
        frame = imgaussfilt(nanstd(stack,[],3),'filtersize',[5,5]);                
    otherwise
        error('unknown var_type');
end

PlotMesoFrame(frame,'caxis',[0 99])

if ~isempty(probe_locations)
    cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','g','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0);
end
end %function