% analyze grab_ne data with auditory cue
clearvars;
close all;
setup_figprop;

root_path = 'V:\HongliWang\grabne_audi\grabne_audi';

logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);
dataIndex = audi_createBehMatFiles(dataIndex);

%% pupil

createPupilFiles(dataIndex);

% plot pupil-tone PSTH

audi_pupil_plots(dataIndex([14],:));
%% fluo
MP_GRAB_preprocess(dataIndex([17,21],:));
audi_GRAB_plots(dataIndex([21],:));

%% check correlation
audi_GRAB_correlation(dataIndex([16,19],:));

audi_pupil_plots(dataIndex([2],:));

% how to find the patchness in the video?
% autocorrealtion?
GRAB_pattern(dataIndex(13,:));
%% plot pupil and fluorescence signal together
% use subject 902

audi_GRAB_pupil_plots(dataIndex(2,:));
audi_GRAB_summary(dataIndex([20,21,23,24,26],:));

%% spontaneous recordings

% spontaneous, GRAB_Ach, before MP. 5,8,10,12
fluo_pupil_plots(dataIndex(1,:));
fluo_pupil_GRAB_summary(dataIndex([1:3],:));

% cross-correlation of pupil and GRAB signals
fluo_pupil_GRAB_correlation(dataIndex([1:17],:));
% Ach
fluo_pupil_corrSummary(dataIndex([1:17],:));

% coherence
coherence_plots(dataIndex([1:3,5,6],:));