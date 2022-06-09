% analyze grab_ne data with auditory cue
clearvars;
close all;
setup_figprop;

root_path = 'Y:\HongliWang\grabne_audi\grabne_audi';

logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);
dataIndex = audi_createBehMatFiles(dataIndex);

%% pupil

createPupilFiles(dataIndex);

% plot pupil-tone PSTH

%audi_pupil_plots(dataIndex([14],:));
%% fluo
MP_GRAB_preprocess(dataIndex([16,25],:));
audi_GRAB_plots(dataIndex([29],:));

%% check correlation
audi_GRAB_correlation(dataIndex([14],:));

% plot pupil PSTH in auditory stimulus trials
audi_pupil_plots(dataIndex([25],:));

%% how to find the patchness in the video?
% autocorrealtion?
GRAB_pattern(dataIndex(13,:));

%% plot pupil and fluorescence signal together
% use subject 902


audi_GRAB_pupil_plots(dataIndex([24:45],:));

savefigpath = fullfile(root_path,'summary');
audi_GRAB_summary(dataIndex,savefigpath);


%% spontaneous recordings

% spontaneous, GRAB_Ach, before MP. 5,8,10,12

% coherence and cross spectrum
fluo_pupil_plots(dataIndex(7,:));
fluo_pupil_GRAB_summary(dataIndex, savefigpath);
% coherence
coherence_plots(dataIndex,savefigpath);

% cross-correlation of pupil and GRAB signals
fluo_pupil_GRAB_correlation(dataIndex([1:18],:));
% Ach
savefigpath = fullfile(root_path,'summary');
fluo_pupil_corrSummary(dataIndex([1:18],:),savefigpath);
