function MP_GRAB_GLM(dataIndex)

%% reference to Sul et al.2011
% running multilinear regression
% model:
% df/f = GLM(choice + outcome + interaction + running_ave + cumulative +
% lick)



nFiles = size(dataIndex,1);


%% go through every session, running GLM separately
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1
        

        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo/LinearRegr');

        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
        
        load(fullfile(fn_fluo.folder,fn_fluo.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        
        
        %saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR.mat']);
        saveRegName = fullfile(savematpath,'glmCR_norm_licknochoice.mat');  % regression for fluo change
        % saveRegName_ITI =  fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);
        
        %% linear regression with C(n+1), observable variables
        
        %  C(n-2) - C(n+1)
        %  R(n-2) - R(n+1)
        %   interaction term
        %if ~exist(saveRegName)
            
            % convert left/right choices into ipsi/contra choices
            % dummy code: ipsi: 0; contra: 1;
            
            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
            
            
            params.xtitle = 'Time from cue (s)';
            params.window = [-3,5];
            params.nback = 0;       %how many trials back to regress against
            params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            params.interaction = false;
            params.ifplot = 0;
            
            params.trigTime(:,1) = trialData.cueTimes - 3;
            params.trigTime(:,2) = trialData.cueTimes + 5;
            %params.trigTime = trialData.cueTimes;
            %only perform analysis on this subset of trials
            
            % 
            choiceTime = trialData.rt+trialData.cueTimes;
            params.eventTime{1} = choiceTime(trialData.response==2);
            params.eventLabel{1} = 'left choice';
            %predictor 2: choice
            params.eventTime{2} = choiceTime(trialData.response==3);
            params.eventLabel{2} = 'right choice';
            params.eventTime{3} = trialData.outcomeTimes;
            params.eventLabel{3} = 'reward';
            
            % remove the first lick from lick times (choice lick)
            leftlickTime = setdiff(sessionData.lickTimes{1,1}, choiceTime);
            rightlickTime = setdiff(sessionData.lickTimes{1,2}, choiceTime);

            params.eventTime{4} = leftlickTime;
            %params.eventTime{4} = sessionData.lickTimes{1,1};
            params.eventLabel{4}= 'left lick';
            params.eventTime{5} = rightlickTime;
            %params.eventTime{5} = sessionData.lickTimes{1,2};
            params.eventLabel{5} = 'right lick';
            tlabel={'left choice','right choice','reward', 'left lick', 'right lick'};
            
            % reg_cr_future=linear_regr( fluo.dia, fluo.t, future_event, params.trigTime, trialMask, params );
            reg_t = cells.t;
            if isfield(cells,'normdFF')
                 reg_dFF = cells.normdFF;
            else
                reg_dFF = cells.dFF;
            end
            hWaitbar = waitbar(0, 'Progress...');
            ifReg = false;
            tic
            %parfor j=1:numel(reg_dFF)
            for j=1:numel(reg_dFF)
                if length(reg_t) > length(reg_dFF{1})
                    [glreg{j}]=gl_regr(reg_dFF{j}, reg_t(1:length(reg_dFF{1})),  params.trigTime, params.eventTime, params, ifReg );
                    %reg_cr{j}=linear_regr( cells.dFF{j}, cells.t(1:length(cells.dFF{1})), future_event, params.trigTime, trialMask, params );
                    
                else
                    %reg_cr{j}=linear_regr( cells.dFF{j}, cells.t, future_event, params.trigTime, trialMask, params );
                    [glreg{j}]=gl_regr(reg_dFF{j}(1:length(reg_t)), reg_t,  params.trigTime, params.eventTime, params,ifReg);
                end
                 progressMessage = sprintf('Progress: %d%%', round(j/numel(reg_dFF)*100));
                 waitbar(j/numel(reg_dFF),hWaitbar, progressMessage);
            end
            toc
            close(hWaitbar);

            all_coeff = [];
            for rr = 1:length(glreg)
                temp_coeff = zeros(1,length(glreg{rr}.t_beta));
                for cc =1:length(glreg{rr}.beta)
                    temp_coeff = [temp_coeff; glreg{1,rr}.beta{1,cc}'];
                end
                all_coeff = cat(3,all_coeff, temp_coeff');
            end
            reg_cr_all.coeff= all_coeff;

            % use bootstrp to get coefficient
            reg_cr_all = getBootstrp(reg_cr_all, 0, 0.05);

            reg_cr_all.regr_time =glreg{1}.t_beta;
            reg_cr_all.numPredictor = length(glreg{1}.beta);
            reg_cr_all.nback = 0;
            reg_cr_all.interaction = false;
            reg_cr_all.pvalThresh= 0.01;

            MP_plot_regrcoef_fluo(reg_cr_all,params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','GLM-norm_lickeffect');    %png format
            saveas(gcf, 'GLM-norm_lickeffect', 'fig');
            saveas(gcf, 'GLM-norm_lickeffect','svg');

            
            all_coeff_nolick = [];
            for rr = 1:length(glreg)
                temp_coeff = zeros(1,length(glreg{rr}.t_beta));
                for cc =1:length(glreg{rr}.beta_reduced)
                    temp_coeff = [temp_coeff; glreg{1,rr}.beta_reduced{1,cc}'];
                end
                all_coeff_nolick = cat(3,all_coeff_nolick, temp_coeff');
            end
            reg_cr_nolick.coeff= all_coeff_nolick;

            % use bootstrp to get coefficient
            reg_cr_nolick = getBootstrp(reg_cr_nolick, 0, 0.05);

            reg_cr_nolick.regr_time =glreg{1}.t_beta;
            reg_cr_nolick.numPredictor = length(glreg{1}.beta_reduced);
            reg_cr_nolick.nback = 0;
            reg_cr_nolick.interaction = false;
            reg_cr_nolick.pvalThresh= 0.01;
            MP_plot_regrcoef_fluo(reg_cr_nolick,params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','GLM-norm_nolick');    %png format
            saveas(gcf, 'GLM-norm_nolick', 'fig');
            saveas(gcf, 'GLM-norm_nolick','svg');
            %%plot the field of view, with significance of choice/reward in every grid
            
            varExplained = [];
            r2_reduced = [];
            for rr = 1:length(glreg)
                r2_reduced = [r2_reduced, glreg{1,rr}.rSquare_reduced];
                varExplained = [varExplained, (glreg{1,rr}.rSquare-glreg{1,rr}.rSquare_reduced)/ glreg{1,rr}.rSquare];
            end

            % plot_regr_fluo(reg_cr,params.pvalThresh,tlabel,params.xtitle);
            % print(gcf,'-dpng','MLR-choiceoutcome_sig_choice0_1');    %png format
            % saveas(gcf, 'MLR-choiceoutcome_sigchoice0_1', 'fig');
            % MP_plot_regrcoef_fluo(reg_cr_future_change,params.pvalThresh,tlabel,params.xtitle);
            
            %         print(gcf,'-dpng','MLR-choiceoutcome_cut_c(n+1)');    %png format
            %         saveas(gcf, 'MLR-choiceoutcome_cut_c(n+1)', 'fig');
            %
 
            
            %% save all regression into one mat file
            % save all these in a structure
            %save(saveRegName, 'reg_cr', 'reg_cr1','reg_cr2','reg_cr3','reg_cr_future','reg_cr_future_ctrl','reg_cr_ctrl');
            save(saveRegName,'glreg', 'glreg_fit');
            %save(saveRegName_ITI,'reg_cr1_change_1','reg_cr2_change_1','reg_cr3_change_1');
            close all;

    end
end
end
