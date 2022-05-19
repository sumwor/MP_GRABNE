% analysis the GRAB_NE signal

clearvars;
close all;
setup_figprop;


%root_path = 'Y:\HongliWang\GRAB_analysis';
%root_path = 'Y:\HongliWang\Ach_analysis';
%root_path = 'C:\Hongli\data\GRAB_analysis';

%root_path = 'V:\HongliWang\NE_analysis';

root_path = 'Y:\HongliWang\NE_analysis_784';

%root_path = 'Y:\HongliWang\Ach_784';

%root_path = 'K:\Ach_784';


%% matching pennies behavior

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('--------------------------------------------------open---------');

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
    MP_GRAB_session(dataIndex.BehPath{ii},dataIndex.LogFileName{ii},savematpath);
end
save_path = fullfile(root_path,'summary','figs_summary');


MP_GRAB_behaviorPerAnimal(dataIndex,save_path);

MP_GRAB_behaviorAll(dataIndex, save_path);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Model fitting

model_path = fullfile(root_path,'mat_models');


MP_GRAB_fittingPerAnimal(dataIndex,model_path);

% save the predicted latent variable into each behavior mat files
MP_GRAB_saveLatent(dataIndex, model_path);

%% fluorescent signal preprocessing
MP_GRAB_preprocess(dataIndex);
% in 910-MP-1213, one img trigger was missing at trial 450 (two trials 450
% + 451) in one imaging file
% solve later

%% simple plots
MP_GRAB_simpleplots(dataIndex);


%% try to make the two baseline equal?
MP_GRAB_checkbaseline(dataIndex);

%% check correlation
audi_GRAB_correlation(dataIndex);
%% fluorescent raw single trial plots
MP_GRAB_rawplots(dataIndex);
%% regression: observable variable, choice, outcome, average reward, cumulative reward
%% save path
model_path = fullfile(root_path,'mat_models');

save_path_fluo = fullfile(root_path,'summary','figs_summary_fluo');

% choice/reward selevitivty

MP_GRAB_selectivity(dataIndex(11,:));
MP_GRAB_prev_selectivity(dataIndex);

% calculate px selectivity
MP_GRAB_selectivity_px(dataIndex);
% make videos of choice/outcome selectivity
MP_GRAB_selectivityVideo(dataIndex);
% spatial pattern
% MP_GRAB_selectivitySpatial(dataIndex);
% MP_GRAB_selectivitySpatial_summary(dataIndex, save_path_fluo);

%% Linear regression with choice and reward
%% for linear regression - load averagve grid intensity, and its corresponding linear regression R2
% running regression (choice and reward) on individual sessions
MP_GRAB_MLR(dataIndex);
% regression results seems unstable, try amount of variance explained?
MP_GRAB_MLR_analysis(dataIndex);
MP_GRAB_MLR_acrossAnimals(dataIndex,save_path_fluo);
MP_GRAB_MLR_acrossSessions(dataIndex([1:9,19:end],:), save_path_fluo)

%%
%
MP_GRAB_MLR_separateSession(dataIndex);
MP_GRAB_MLR_separateSummary(dataIndex, save_path_fluo);

% latent variable
MP_GRABRL_MLR(dataIndex([10:17],:))
MP_GRABRL_MLR_acrossSessions(dataIndex([1:9,19:end],:), save_path_fluo)

MP_GRABRL_sumQ_MLR(dataIndex[10:17],:))
MP_GRABRL_sumQ_MLR_acrossSessions(dataIndex, save_path_fluo)

MP_GRABRL_RPE_MLR(dataIndex);
MP_GRABRL_RPE_MLR_acrossSessions(dataIndex([1:9,19:end],:), save_path_fluo)

MP_GRAB_selectivitySpatial(dataIndex);
MP_GRAB_selectivityCorr(dataIndex);
MP_GRAB_selectivitySummary(dataIndex([1:9,18:end],:), save_path_fluo);

MP_GRAB_temporalCorr(dataIndex);

% temporal correlation (same grid, correlation of different variables)


%% random forest is better in estimating weak signals
MP_GRAB_RF(dataIndex);
MP_GRAB_RF_acrossSessions(dataIndex, save_path_fluo)

MP_GRAB_RF_RL(dataIndex);
MP_GRAB_RF_RL_acrossSessions(dataIndex, save_path_fluo)

% motion kernel?
% df/f - reward rate/latent variable (sum of Q etc. )
MP_GRAB_tonic(dataIndex)

%% GLM
MP_GRAB_GLM_findlambda(dataIndex);

lambda = 0.2;  % manually picking lambda = 0.2 for every session for now
MP_GRAB_GLM(dataIndex,lambda);
% bayesian linear regression
MP_GRAB_BLM(dataIndex);
%% PCA
MP_GRAB_PCA(dataIndex);

%% clustering
MP_GRAB_clustering(dataIndex);

%% find the patchness?
