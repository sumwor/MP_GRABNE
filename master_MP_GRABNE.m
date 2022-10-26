% analysis the GRAB_NE signal

clearvars;
close all;
setup_figprop;


root_path_NE = 'V:\HongliWang\GRAB\NE_analysis_784';
root_path_ACh = 'V:\HongliWang\GRAB\ACh_analysis_784';
%root_path = 'K:\Ach_784';


%% matching pennies behavior

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('--------------------------------------------------open---------');

% Look for data files and create a database index

%% NE
logfilepath_NE = fullfile(root_path_NE,'data');
analysispath_NE = fullfile(root_path_NE,'analysis');
dataIndex_NE = makeDataIndex(logfilepath_NE, analysispath_NE);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex_NE = MP_GRAB_createBehMatFiles(dataIndex_NE,'NE');

% sort Index according to experiment data

dataIndex_NE = sortdataIndex(dataIndex_NE);

% Determine if each session fulfill performance criteria (Matching Pennies)
MP_determineBehCriteria(dataIndex_NE([1:13],:));

% ACh
logfilepath_ACh = fullfile(root_path_ACh,'data');
analysispath_ACh = fullfile(root_path_ACh,'analysis');
dataIndex_ACh = makeDataIndex(logfilepath_ACh, analysispath_ACh);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex_ACh = MP_GRAB_createBehMatFiles(dataIndex_ACh,'ACh');

% sort Index according to experiment data

dataIndex_ACh = sortdataIndex(dataIndex_ACh);

% Determine if each session fulfill performance criteria (Matching Pennies)
MP_determineBehCriteria(dataIndex_ACh);

%% behavior

% NE
nFiles_NE = size(dataIndex_NE([1:13],:),1);

for ii = 1:nFiles_NE
    savematpath_NE = dataIndex_NE.BehPath{ii};
    MP_GRAB_session(dataIndex_NE.BehPath{ii},dataIndex_NE.LogFileName{ii},savematpath_NE);
end
save_path_NE = fullfile(root_path_NE,'summary','figs_summary');


MP_GRAB_behaviorPerAnimal(dataIndex_NE,save_path_NE);

MP_GRAB_behaviorAll(dataIndex_NE, save_path_NE);

% ACh
nFiles_ACh = size(dataIndex_ACh,1);

for ii = 1:nFiles_ACh
    savematpath_ACh = dataIndex_ACh.BehPath{ii};
    MP_GRAB_session(dataIndex_ACh.BehPath{ii},dataIndex_ACh.LogFileName{ii},savematpath_ACh);
end
save_path_ACh = fullfile(root_path_ACh,'summary','figs_summary');


MP_GRAB_behaviorPerAnimal(dataIndex_ACh,save_path_ACh);

MP_GRAB_behaviorAll(dataIndex_ACh, save_path_ACh);

% summary
MP_GRAB_behComp(dataIndex_ACh, dataIndex_NE, save_path_ACh)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Model fitting

% NE
model_path_NE = fullfile(root_path_NE,'mat_models');


MP_GRAB_fittingPerAnimal(dataIndex_NE,model_path_NE);

% save the predicted latent variable into each behavior mat files
MP_GRAB_saveLatent(dataIndex_NE, model_path_NE);

%ACh

model_path_ACh = fullfile(root_path_ACh,'mat_models');


MP_GRAB_fittingPerAnimal(dataIndex_ACh,model_path_ACh);

% save the predicted latent variable into each behavior mat files
MP_GRAB_saveLatent(dataIndex_ACh, model_path_ACh);
%% fluorescent signal preprocessing

% NE



%ACh
MP_GRAB_preprocess(dataIndex_ACh);
% in 910-MP-1213, one img trigger was missing at trial 450 (two trials 450
% + 451) in one imaging file
% solve later

%% simple plots
MP_GRAB_simpleplots(dataIndex_NE);
MP_GRAB_simpleplots(dataIndex_ACh);

%% try to make the two baseline equal?
MP_GRAB_checkbaseline(dataIndex_NE);
MP_GRAB_checkbaseline(dataIndex_ACh);
%% check correlation
audi_GRAB_correlation(dataIndex);
%% fluorescent raw single trial plots
MP_GRAB_rawplots(dataIndex);

%% regression: observable variable, choice, outcome, average reward, cumulative reward
%% save path

% NE
model_path_NE = fullfile(root_path_NE,'mat_models');

save_path_fluo_NE = fullfile(root_path_NE,'summary','figs_summary_fluo');

% ACh
model_path_ACh = fullfile(root_path_ACh,'mat_models');

save_path_fluo_ACh = fullfile(root_path_ACh,'summary','figs_summary_fluo');
%

%
% % choice/reward selevitivty
%
% MP_GRAB_selectivity(dataIndex(11,:));
% MP_GRAB_prev_selectivity(dataIndex);
%
% % calculate px selectivity
% MP_GRAB_selectivity_px(dataIndex);
% % make videos of choice/outcome selectivity
% MP_GRAB_selectivityVideo(dataIndex);
% % spatial pattern
% % MP_GRAB_selectivitySpatial(dataIndex);
% % MP_GRAB_selectivitySpatial_summary(dataIndex, save_path_fluo);

%% Linear regression with choice and reward
%% for linear regression - load averagve grid intensity, and its corresponding linear regression R2
% running regression (choice and reward) on individual sessions

% NE

MP_GRAB_MLR(dataIndex_NE([1:9],:));
% regression results seems unstable, try amount of variance explained?
MP_GRAB_MLR_analysis(dataIndex_NE);
MP_GRAB_MLR_acrossAnimals(dataIndex_NE,save_path_fluo_NE);
MP_GRAB_MLR_acrossSessions(dataIndex_NE, save_path_fluo_NE)

MP_GRAB_MLR_separateSession(dataIndex_NE);
MP_GRAB_MLR_separateSummary(dataIndex_NE, save_path_fluo_NE);

% ACh
MP_GRAB_MLR(dataIndex_ACh);

% bootstrap not useful, estimation still unsmooth
% MP_GRAB_MLR_bootstrap(dataIndex_ACh);
% regression results seems unstable, try amount of variance explained?
MP_GRAB_MLR_analysis(dataIndex_ACh);
MP_GRAB_MLR_acrossAnimals(dataIndex_ACh,save_path_fluo_ACh);
MP_GRAB_MLR_acrossSessions(dataIndex_ACh, save_path_fluo_ACh);

MP_GRAB_MLR_separateSession(dataIndex_ACh);
MP_GRAB_MLR_separateSummary(dataIndex_ACh, save_path_fluo_ACh);


% compare ACh and NE outcome curves
MP_GRAB_MLR_comp(dataIndex_ACh,dataIndex_NE([1:19,28:end],:));

%% latent variable

%NE
MP_GRABRL_MLR(dataIndex_NE([1:9],:))
MP_GRABRL_MLR_acrossSessions(dataIndex_NE, save_path_fluo_NE)

MP_GRABRL_sumQ_MLR(dataIndex_NE([10:17],:))
MP_GRABRL_sumQ_MLR_acrossSessions(dataIndex_NE, save_path_fluo_NE)

MP_GRABRL_RPE_MLR(dataIndex_NE([1:9],:));
MP_GRAB_RPE_MLR_acrossAnimals(dataIndex_NE, save_path_fluo_NE)
MP_GRABRL_RPE_MLR_acrossSessions(dataIndex_NE, save_path_fluo_NE)

% ACh
MP_GRABRL_MLR(dataIndex_ACh)
MP_GRABRL_MLR_acrossSessions(dataIndex_ACh, save_path_fluo_ACh)

MP_GRABRL_sumQ_MLR(dataIndex_ACh)
MP_GRABRL_sumQ_MLR_acrossSessions(dataIndex_ACh, save_path_fluo_ACh)

MP_GRABRL_RPE_MLR(dataIndex_ACh);
MP_GRAB_RPE_MLR_acrossAnimals(dataIndex_ACh, save_path_fluo_ACh)
MP_GRABRL_RPE_MLR_acrossSessions(dataIndex_ACh, save_path_fluo_ACh)

%% spatial & temporal patterns

% NE
MP_GRAB_selectivitySpatial(dataIndex_NE);
MP_GRAB_selectivityCorr(dataIndex_NE);
MP_GRAB_selectivitySummary(dataIndex_NE, save_path_fluo_NE);

MP_GRAB_clusterEval(dataIndex_NE([1:9,18:end],:), save_path_fluo_NE, save_path_mat_NE)

MP_GRAB_temporalCorr(dataIndex_NE([1:9],:));
save_path_mat_NE = fullfile(root_path_NE,'summary','data_summary');
MP_GRAB_temporalCorrSummary(dataIndex_NE,save_path_fluo_NE,save_path_mat_NE)
% temporal correlation (same grid, correlation of different variables)

% ACh
MP_GRAB_selectivitySpatial(dataIndex_ACh);
MP_GRAB_selectivityCorr(dataIndex_ACh);
MP_GRAB_selectivitySummary(dataIndex_ACh, save_path_fluo_ACh);

MP_GRAB_clusterEval(dataIndex_ACh, save_path_fluo_ACh, save_path_mat_ACh)

MP_GRAB_temporalCorr(dataIndex_ACh);
save_path_mat_ACh = fullfile(root_path_ACh,'summary','data_summary');
MP_GRAB_temporalCorrSummary(dataIndex_ACh,save_path_fluo_ACh, save_path_mat_ACh)


% compare NE and ACh result
MP_GRAB_tempComp(save_path_mat_NE, save_path_mat_ACh,save_path_fluo_ACh)

%% random forest is better in estimating weak signals
% MP_GRAB_RF(dataIndex);
% MP_GRAB_RF_acrossSessions(dataIndex, save_path_fluo)
%
% MP_GRAB_RF_RL(dataIndex);
% MP_GRAB_RF_RL_acrossSessions(dataIndex, save_path_fluo)
%
% % motion kernel?
% % df/f - reward rate/latent variable (sum of Q etc. )
% MP_GRAB_tonic(dataIndex)
%
% %% GLM
% MP_GRAB_GLM_findlambda(dataIndex);
%
% lambda = 0.2;  % manually picking lambda = 0.2 for every session for now
% MP_GRAB_GLM(dataIndex,lambda);
% % bayesian linear regression
% MP_GRAB_BLM(dataIndex);
%% PCA
MP_GRAB_PCA(dataIndex);

%% clustering
MP_GRAB_clustering(dataIndex);

%% find the patchness?
