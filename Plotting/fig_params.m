classdef fig_params  
    properties
        %Line Plot Options
        p_line_width = 1.5; 
        
        %distribution line plot options
        dl_line_width = 1.5; 
        dl_alpha = 0.5;
        
        %violin plots
        vp_alpha = 0.3;
        vp_dist_w = 0.4;
        
        %colors
        c_discovery = [0.3 0.7 1]; %light blue
        c_pca = [0.14 0.6412 0.2471] %dark green 0.1412    0.5412    0.2471
        c_nmf = [0.850,0.3723,0.00784]; %orange
        c_wh_same = [0 0 0.75]; %withheld, same animal
        c_wh_diff = [0.3 0.7 .4]; %withheld, different animal
        c_static = [0.5 0.5 0.5]; %static 
        c_sensory = [0.3 0.7 1; 0.7,0.1,0; 0.3 0.7 0] %air,vis,aud. 
        c_visual = [0, 0, 1; 0.7,0.1,0 ];
        c_tactile = [0, 0, 0.75; 0.3, 0.7, 1]; %[1, 0, 0; 0.3, 0.7, 1];
        c_spacetime = [0.7608 0 0.0627]; %the space, time decomposition
        
        %Global figure options
        font_size = 16;
        font_name = 'Arial';
        font_weight = 'normal';
        units = 'centimeters';
        line_width = 1;            
        default_color = [0.5 0.5 0.5]; %for any default color plots. 
        
    end

    methods
        function FormatAxes(obj,ax_handle)
            set(ax_handle,'fontsize',obj.font_size,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'linewidth',obj.line_width,...
                'fontsize',obj.font_size,...
                'box','off');
        end %end function
        
        function SetTitle(obj,ax_handle,title_str)
            title(ax_handle,title_str,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'fontsize',obj.font_size)
        end%end function
        
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











