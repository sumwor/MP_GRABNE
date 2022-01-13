% analyze grab_ne data with auditory cue
clearvars;
close all;
setup_figprop;

root_path = 'E:\grabne_audi\grabne_audi';

logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);
dataIndex = audi_createBehMatFiles(dataIndex);

%% pupil

createPupilFiles(dataIndex);

% plot pupil-tone PSTH
audi_pupil_plots(dataIndex([14],:));
%% fluo
MP_GRAB_preprocess(dataIndex([10,14],:));
audi_GRAB_plots(dataIndex([10,14],:));

%% check correlation
audi_GRAB_correlation(dataIndex([10,14],:));

% how to find the patchness in the video?
% autocorrealtion?

%% plot pupil and fluorescence signal together
% use subject 902
audi_GRAB_pupil_plots(dataIndex([14],:));
%% spontaneous recordings
fluo_pupil_plots(dataIndex(10,:));
% coherence
coherence_plots(dataIndex([1:3,5,6],:));