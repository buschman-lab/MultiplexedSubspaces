function BatchBehavioralProcessing()

[fn,~] = GrabFiles('\w*_parsedVideo.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BehavioralData\'});

if numel(fn)<6
    fn = cell(1,6);
    for cur_rec = 1:6
       fn{cur_rec}=ParseBehavioralVideo(cur_rec);
    end
end

%could also probably run on a parpool
for cur_rec = 1:6
    BehavProcessing(cur_rec,fn{cur_rec});
end

fprintf('\nDONE WITH EVERYTHING')

end %function end

% do locally, spock sucks for video. Hangs all the time (probably could
% resolve with different codecs... but that's for another day). 

% %send off to spock for processing
% %Camden MacDowell
% parameter_class = 'general_params_corticaldynamics';
% % Open ssh connection
% username = input(' Spock Username: ', 's');
% password = passcode();
% s_conn = ssh2_config('spock.princeton.edu',username,password);
% 
% chunk = load(fn{1},'chunk');
% chunk = chunk.chunk;
% for cur_rec = 1:6
%     for cur_c = 1%:2%numel(chunk)
%     script_name = WriteBashScript(parameter_class,sprintf('%d',1),'BehavProcessing',{cur_c,cur_rec,ConvertToBucketPath(fn{cur_rec})},...
%         {'%d','%d',"'%s'"},...
%         'sbatch_time',5,'sbatch_memory',24,...
%         'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/BehavioralAnalysisFunctions/");
% 
%     ssh2_command(s_conn,...
%     ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
%     sprintf('sbatch %s',script_name)]);  
%     end
% end
% 
% %close out connection
% ssh2_close(s_conn);
% clear username password s_conn
