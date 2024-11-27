function MP_GRAB_MLR_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
setup_figprop
set(0, 'DefaultFigureRenderer', 'painters');
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

% load the first file to get some parameters
all_coeff_future = [];
all_pval_future = [];

all_coeff_future_nolick = [];
all_pval_future_nolick = [];

percentage_outcome_session = [];
percentage_outcome_session_nolick = [];
percentage_choice_session = [];
percentage_choice_session_nolick = [];

pvalThresh = 0.01;
for ii = 1:nFiles
    display(ii);
    temp_all_pval_future = [];
    
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
    % load behavior files
       fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    
     %saveRegName = fullfile(savematpath,'regCR_norm.mat');  % regression for fluo change
    saveRegName = fullfile(savematpath,'regCR_norm_lick.mat');
    if exist(saveRegName)
        load(saveRegName)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
    
    % load the MLR with C(n+1)
        for rr = 1:length(reg_cr) 
            all_coeff_future = cat(3,all_coeff_future, reg_cr{rr}.coeff);
            all_pval_future = cat(3, all_pval_future, reg_cr{rr}.pval);
            temp_all_pval_future = cat(3, temp_all_pval_future, reg_cr{rr}.pval);
        end
        reg_all.regr_time = reg_cr{1}.regr_time;
        reg_all.numPredictor = reg_cr{1}.numPredictor;
        reg_all.nback = reg_cr{1}.nback;
        reg_all.interaction = reg_cr{1}.interaction;
        
        nCells = length(reg_cr);
        temp_percent_outcome = 100*sum(temp_all_pval_future(:,7,:)<pvalThresh,3)/nCells;
        percentage_outcome_session = cat(3, percentage_outcome_session,  temp_percent_outcome);
        temp_percent_choice = 100*sum(temp_all_pval_future(:,3,:)<pvalThresh,3)/nCells;
        percentage_choice_session = cat(3, percentage_choice_session,  temp_percent_choice);

    end
end
        

% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

%% load nolick regression result
for ii = 1:nFiles
    display(ii);
    temp_all_pval_future = [];
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
    % load behavior files
       fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    
     %saveRegName = fullfile(savematpath,'regCR_norm.mat');  % regression for fluo change
    saveRegName = fullfile(savematpath,'regCR_norm.mat');
    if exist(saveRegName)
        load(saveRegName)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
          
    % load the MLR with C(n+1)
        for rr = 1:length(reg_cr) 
            all_coeff_future_nolick = cat(3,all_coeff_future_nolick, reg_cr{rr}.coeff);
            all_pval_future_nolick = cat(3, all_pval_future_nolick, reg_cr{rr}.pval);
            temp_all_pval_future = cat(3, temp_all_pval_future, reg_cr{rr}.pval);
        end
        reg_all.regr_time = reg_cr{1}.regr_time;
        reg_all.numPredictor = reg_cr{1}.numPredictor;
        reg_all.nback = reg_cr{1}.nback;
        reg_all.interaction = reg_cr{1}.interaction;
        
        nCells = length(reg_cr);
        temp_percent_outcome = 100*sum(temp_all_pval_future(:,7,:)<pvalThresh,3)/nCells;
        percentage_outcome_session_nolick = cat(3, percentage_outcome_session_nolick,  temp_percent_outcome);
        temp_percent_choice_nolick = 100*sum(temp_all_pval_future(:,3,:)<pvalThresh,3)/nCells;
        percentage_choice_session_nolick = cat(3, percentage_choice_session_nolick,  temp_percent_choice_nolick);

     % load the ITI regression (n+1 and n)

    end
end


%% stats between lick and nolick

percentage_outcome_session = squeeze(percentage_outcome_session);
percentage_outcome_session_nolick = squeeze(percentage_outcome_session_nolick);
percentage_choice_session = squeeze(percentage_choice_session);
percentage_choice_session_nolick = squeeze(percentage_choice_session_nolick);

[p, table] = anova_rm({percentage_outcome_session' percentage_outcome_session_nolick'});
[p, table] = anova_rm({percentage_choice_session' percentage_choice_session_nolick'});

%%  linear regression with C(n+1)
reg_cr_all.coeff= all_coeff_future;

% use bootstrp to get coefficient
reg_cr_all = getBootstrp(reg_cr_all, 0, 0.05);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);
reg_cr_all.regr_time = reg_all.regr_time;
reg_cr_all.numPredictor = reg_all.numPredictor;
reg_cr_all.nback = reg_all.nback;
reg_cr_all.interaction = reg_all.interaction;
reg_cr_all.pvalThresh= 0.01;
% 
% xtickangle(45)
% ylabel('Coefficients (a.u.)');
% title('Coefficient for pupil change - choice and reward');
xtitle='Time from cue (s)';
 tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
        'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward', 'contralick', 'ipsilick'};
pvalThresh=0.01;
MP_plot_regrcoef_fluo(reg_cr_all,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-norm-choiceoutcome_future_lick_averageSession');    %png format
saveas(gcf, 'MLR-norm-choiceoutcome_future_lick_averageSession', 'fig');
saveas(gcf, 'MLR-norm-choiceoutcome_future_lick_averageSession','svg');

% plot the figure as number of session that is significant
reg_sig.coeff = all_coeff_future;
reg_sig.pval = all_pval_future;
reg_sig.regr_time = reg_all.regr_time;
reg_sig.numPredictor = reg_all.numPredictor;
reg_sig.nback = reg_all.nback;
reg_sig.interaction = reg_all.interaction;
reg_sig.pvalThresh= 0.01;

reg_sig_nolick.coeff = all_coeff_future_nolick;
reg_sig_nolick.pval = all_pval_future_nolick;

% test for outcome
nCells = size(reg_sig.coeff,3);
pvalThresh = 0.01;
percentage_outcome = 100*sum(reg_sig.pval(:,7,:)<pvalThresh,3)/nCells;
percentage_outcome_nolick = 100*sum(reg_sig_nolick.pval(:,7,:)<pvalThresh,3)/nCells;

percentage_choice = 100*sum(reg_sig.pval(:,3,:)<pvalThresh,3)/nCells;
percentage_choice_nolick = 100*sum(reg_sig_nolick.pval(:,3,:)<pvalThresh,3)/nCells;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);

if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr_fluo(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr_fluo(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-norm-choiceoutcome_future_lick_sigSession');    %png format
saveas(gcf, 'MLR-norm-choiceoutcome_future_lick_sigSession', 'fig');
saveas(gcf, 'MLR-norm-choiceoutcome_future_lick_sigSession','svg');

% save the analysis
savematpath = fullfile(savefigpath,'Result');
if ~exist(savematpath)
    mkdir(savematpath);
end
save(fullfile(savematpath,'regression1Sum'),'reg_cr_all','reg_sig');
%%
close all

end