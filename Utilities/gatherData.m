function params = gatherData(numIter,base,grpflag,directory,variables)
    fprintf('\n Compiling data...\n');
    %Go to the directory and load the recordings, divide up by group and    
    cd(directory)
    
    load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO\MouseNumInfoByRec.mat')
    %Get grp info
    grp = isVPA(mousenum);
    if grpflag == 1 %do vpa animals
        rmvFlag = 0; %remove Saline
    elseif grpflag ==0 %do saline animals
        rmvFlag = 1; %remove VPA
    else %Keep both
        rmvFlag =2; 
    end
    
    %Load variables of interest
    params = struct();
    COUNT =1;
    for i = 1:numIter
        if grp(i)~=rmvFlag
            try
                temp = load(sprintf('%s_%d.mat',base,i),variables{:});
                for j= 1:numel(variables)
                    params(COUNT).(variables{j}) = temp.(variables{j});
                end            
                %Add mouse number
                params(COUNT).Mouse = mousenum(i);
            catch
                for j= 1:numel(variables)
                    params(COUNT).(variables{j}) = NaN;                
                end 
            end
            COUNT = COUNT+1;
        end
    end
end