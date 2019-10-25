classdef fig_params  
    properties
        %Line Plot Options
        p_line_width = 1.5; 
        p_smooth_window = 30; %best if set to block size
        
        
        %Scatter Plot Options
        point_color = [0.75 0.75 0.75];
        point_size = 10; 
        point_marker_repitition = 'o';
        point_repition_size = 12;
        
        
        %patch plots
        patch_alpha = 0.25;
        patch_alpha_min = 0.25;
        patch_alpha_max = 0.25;
        
        %quiver plot 
        q_window = 5; %number of samples to smooth for plotting shift 
        q_head_length = 3;
        q_head_width = 3;
        q_scale = 0.9 %scale (0 to 1) of the length of the line between points
        
        %Target Specific Colors
        target_colors = [0.98 0.41 0.17;0.3 0.2 1];
        target_marker_scale = 10; %how much larger make the target marker than the others
        target_learn_colors = [0.75 0.1 0;0 0.1 0.75];        
        
        %bar plots
        face_alpha = 0.75;
        edge_alpha = 1;
        
        %Global figure options
        font_size = 16;
        font_name = 'Arial';
        font_weight = 'normal';
        units = 'centimeters';
        line_width = 1.5;            
        
        
    end

    methods
        function FormatAxes(obj,ax_handle)
            set(ax_handle,'fontsize',obj.font_size,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'linewidth',obj.line_width,...
                'box','off');
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
            end
            function SaveFigs(obj,handles,filetype,name,savedir,fancyflag)                
                saveCurFigs(handles,filetype,name,savedir,fancyflag)
                close(handles)
            end
    end

end











