function [stats,group] = CompileStats(in_fns,variables,group_ids,group)
    if nargin <4; group =[]; end %set to 'none' for no group selection
    
    if strcmp(group,'none')==1
        %don't do anything
    elseif ~isempty(group) %split by predefined groups according to group_ids
        in_fns(~ismember(group,group_ids))=[];
    else %split by the predefined groups
        fprintf('\n Loading mouse numbers. This can take a while... ')
        tic
        mouse_num = cellfun(@(x) load(x,'MouseNum'),in_fns,'UniformOutput',0);        
        mouse_num = cellfun(@(x) str2num(cell2mat(regexp(x.MouseNum,'\d+','match'))),mouse_num,'UniformOutput',0);
        toc
        group = isVPA(cat(2,mouse_num{:}));
        in_fns(~ismember(group,group_ids))=[];        
        fprintf('\n Done loadings mouse numbers')
    end 
    
    fprintf('\n Compiling data...\n');
    %Load variables of interest
    stats = struct();
    for i = 1:numel(in_fns)
%         try
        if mod(i,round(0.1*numel(in_fns))) ==0
            fprintf('\t%g%% Complete\n', round(i./numel(in_fns)*100,2));
        end
        temp = load(in_fns{i},variables{:});
        for j= 1:numel(variables)
            stats(i).(variables{j}) = temp.(variables{j});
        end  
%         catch
%         end
    end
        
end %function end 


% scratch
%         mouse_num = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO\MouseNumInfoByRec.mat');
%         group = isVPA(mouse_num.mousenum);