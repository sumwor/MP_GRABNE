% analyze grab_ne data with auditory cue
clearvars;
close all;
setup_figprop;

root_path = 'E:\data\grabne_audi';

logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);
dataIndex = audi_createBehMatFiles(dataIndex);

%% pupil

createPupilFiles(dataIndex);

% plot pupil-tone PSTH
audi_pupil_plots(dataIndex);
%% fluo
MP_GRAB_preprocess(dataIndex);
audi_GRAB_plots(dataIndex);

