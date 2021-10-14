% analysis the GRAB_NE signal 

clearvars;
close all;
setup_figprop;

root_path = 'Y:\HongliWang\GRAB_analysis';
%% matching pennies behavior 

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('-----------------------------------------------------------');

% Look for data files and create a database index
logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex = MP_createBehMatFiles(dataIndex);

% sort Index according to experiment data

dataIndex = sortdataIndex(dataIndex);

% Determine if each session fulfill performance criteria (Matching Pennies)
MP_determineBehCriteria(dataIndex);

%% behavior
nFiles = size(dataIndex,1);

for ii = 1:nFiles
    savematpath = dataIndex.BehPath{ii};
    MP_session(dataIndex.BehPath{ii},dataIndex.LogFileName{ii},savematpath);
end

%% fluorescent signal preprocessing
MP_GRAB_preprocess(dataIndex);

%% simple plots
MP_GRAB_simpleplots(dataIndex);

%% model fitting

%% fluorescent raw single trial plots
MP_GRAB_rawplots(dataIndex);
%% regression: observable variable, choice, outcome, average reward, cumulative reward
%% save path
model_path = fullfile(root_path,'mat_models');

save_path_fluo = fullfile(root_path,'summary','figs_summary_fluo');

% choice/reward selevitivty
MP_GRAB_selectivity(dataIndex);
%% Linear regression with choice and reward
% running regression (choice and reward) on individual sessions
MP_GRAB_MLR(dataIndex);
% regression results seems unstable, try amount of variance explained?
MP_GRAB_MLR_analysis(dataIndex);
MP_GRAB_MLR_acrossSessions(dataIndex, save_path_fluo)

%% random forest is better in estimating weak signals
MP_GRAB_RF(dataIndex);
MP_GRAB_RF_acrossSessions(dataIndex, save_path_fluo)

%% GLM
MP_GRAB_GLM(dataIndex);
% bayesian linear regression
MP_GRAB_BLM(dataIndex);
%% PCA
MP_GRAB_PCA(dataIndex);