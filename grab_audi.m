% analyze grab_ne data with auditory cue
clearvars;
close all;
setup_figprop;

root_path = 'E:\data\grabne_audi\grabne_audi';

logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);
dataIndex = audi_createBehMatFiles(dataIndex);

%% pupil

createPupilFiles(dataIndex);

% plot pupil-tone PSTH
audi_pupil_plots(dataIndex(4:6,:));
%% fluo
MP_GRAB_preprocess(dataIndex);
audi_GRAB_plots(dataIndex(11,:));

%% check correlation
audi_GRAB_correlation(dataIndex);

% how to find the patchness in the video?
% autocorrealtion?

%% plot pupil and fluorescence signal together
% use subject 902
audi_GRAB_pupil_plots(dataIndex(8,:));
%% spontaneous recordings
fluo_pupil_plots(dataIndex([1:3,5,6],:));
% coherence
coherence_plots(dataIndex([1:3,5,6],:));