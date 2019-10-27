function script_name = WriteBashScript(uniqID,func_name,input_val,input_type)
%Writes a spock bash script that will run func_name with variables defined
%by input_val. 
%input_val and input_type are equal length cell arrays. 
%example: input_type = {'"%s"'}; sprintf(input_type{1},'$SLURM_ARRAY_TASK_ID'); 

gp = general_params;
%Get the text from the base
text = fileread('spock_base.sh');

%keep this local for right now;
script_name = sprintf('run_%s.sh',uniqID);
file_path = [gp.local_bucket gp.dynamic_script_path];
if ~exist(file_path)
    mkdir(file_path);
end
fid = fopen([file_path script_name], 'wt');

fprintf(fid,'%s','#!/usr/bin/env bash'); %add the first line of the header

%give the script a funny name
temp = EntertainingSpockNames;
fprintf(fid,"\n%s'%s'",'#SBATCH -J ',temp);

%Add the base script
fprintf(fid,'%s',text);

%Add the specific function call
try %make sure to close the fid even if crash
    %Create the variable lengthed inputs
    temp = {};
    for i = 1:numel(input_val)
        if i~=numel(input_val)
            temp{i} = [sprintf([input_type{i}],input_val{i}),','];
        else %final round, no comma
            temp{i} = sprintf([input_type{i}],input_val{i});
        end
    end   
    fprintf(fid, ['xvfb-run -d matlab -nosplash -nodisplay -nodesktop -r ',...
        sprintf('"try;%s(',func_name),... %the function
        sprintf('%s',[temp{:}]),...
        ');catch me;disp(me.message);end;exit;"']); %the inputs              
    fclose(fid);
    
    %convert to unix
    unix2dos([file_path script_name],1)
catch     
    fclose(fid);
    error('Failed generating bash script'); 
end

