function Spock_ClassifyBehavioralWaveForms(filename)
rng('default') %for reproducibility
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end
disp(1)
cd('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/Spock_Temp_Data/');

load(filename);
[~, Observed, ~] = SVMClassifier_Binary(savedata,'featureselect','none',...
'nshuf',0,'kernel','linear','optimize',1,'verbose',1);

save(['Processed_' filename],'Observed');

end