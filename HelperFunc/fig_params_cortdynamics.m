classdef fig_params_cortdynamics  
    properties
        %Line Plot Options
        p_line_width = 1.5; 
        
        %distribution line plot options
        dl_line_width = 1.5; 
        dl_alpha = 0.5;
        
        %violin plots
        vp_alpha = 0.5;
        vp_dist_w = 0.4;
        
        %markers
        markersizesmall = 5; 
        markersizebig = 7; 
        markers = {'o','d','*','x','v','+','^','s'};


        %colors
        c_area = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];...
                [0.4 0.76 0.9]; [.4 .2 .7]; [.7 .2 .1]];
        c_none = [0.5 0.5 0.5]; %light blue
        c_lr = [0.3 0.7 1]; %light blue
        c_glm = [0.14 0.6412 0.2471] %dark green 0.1412    0.5412    0.2471
        c_ff = [0.850,0.3723,0.00784]; %orange
        
        %Global figure options
        font_size = 12;
        font_name = 'Arial';
        font_weight = 'normal';
        units = 'centimeters';
        line_width = 1;            
        default_color = [0.5 0.5 0.5]; %for any default color plots. 
        
        %other
        noisemotif = 2; %the one motif that corresponds to noise
        
    end

    methods
        function FormatAxes(obj,ax_handle)
            set(ax_handle,'fontsize',obj.font_size,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'linewidth',obj.line_width,...
                'fontsize',obj.font_size,...
                'box','off','GridAlpha',0.1);
        end %end function
        
        function SetTitle(obj,ax_handle,title_str)
            title(ax_handle,title_str,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'fontsize',obj.font_size)
        end%end function
        
        function FormatLSline(obj,ax_handle)
            delete(ax_handle(1))
            for i = 2:4
                ax_handle(i).LineWidth=2;
                ax_handle(i).Color=[0.5 0.5 0.5];
            end
            legend off
        end       
        
        function FigureSizing(obj,handles,ax_position,fig_position)
            for i = 1:numel(handles)
                set(0, 'currentfigure', handles(i)); 
                set(gca,'units','centimeters',...
                    'position',ax_position)
                if ~isempty(fig_position)
                    set(handles(i),'units','centimeters',...
                    'position',fig_position)
                end
            end
        end%end function
        
        function SaveFigs(obj,handles,filetype,name,savedir,fancyflag)                
            saveCurFigs(handles,filetype,name,savedir,fancyflag)
            close(handles)
        end%end function
        
        function col = GenDefaultColorArray(obj,array_size)
            col = cell(1,array_size);
            for i = 1:array_size
                col{i} = obj.default_color;
            end
        end %end function
       
       
        
    end

end











