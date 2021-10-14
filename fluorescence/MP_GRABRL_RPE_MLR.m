function MP_GRABRL_RPE_MLR(dataIndex)

% add the cutpoint later

nFiles = size(dataIndex,1);
%% reference: Sul et al.
% s(t) = a0 + a1C(t) + a2R(t) + a3X(t) + a4deltaQ(t) + a5Qc(t) + A(t)
% no A(t) first
%% fit the model per session
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_FQRPECKlatentV.mat']);
    load(fn_latent);
    % load pupil files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    
    if length(fn_fluo) == 1
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
        
        load(fullfile(fn_fluo.folder,fn_fluo.name));
        
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        %saveMLRmatpath = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_CK.mat']);
        saveRegName = fullfile(savematpath,'regRL_RPE.mat');  % regression for fluo changehange
        % saveMLRmatpath_outcome = fullfile(dataIndex.BehPath{ii},[fn_beh.name(1:end-7),'regRL_lag0_outcome_cut_fitall.mat']);
        
        if ~exist(saveRegName)
            params=[];
            
           choice = NaN(size(trials.left));
        choice(trials.left) = 0;
        choice(trials.right) = 1;
        % dummycode left: 0, right 1
        %C(n)
        params.trigEvent = choice;
        
        
        params.trigEvent2= [NaN; params.trigEvent(1:end-1,1)]; % c(n-1)
         
        %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
        reward = NaN(size(trials.go));
        reward(trials.reward) = 1;
        reward(trials.noreward) = 0;
        params.trigEvent3= [NaN;reward(1:end-1)];   % r(n-1)
        
        % dQ
        params.trigEvent4 = stats_new.ql-stats_new.qr; % delta q
        %RPE
        params.trigEvent5= stats_new.rpe;
        % dK
        params.trigEvent6 = stats_new.ckl-stats_new.ckr;
        % CKE
        % 1-choice choice
        params.trigEvent7 = NaN(size(trials.go));
        params.trigEvent7(stats_new.c==-1) = 1-stats_new.ckl(stats_new.c==-1); % choosing left
        params.trigEvent7(stats_new.c==1) = 1-stats_new.ckr(stats_new.c==1);
        
        params.trigEvent8 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent8(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent8(kk) = sum(trials.reward(kk-19:kk))/20;
            end
        end
        
        params.trigEvent9=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent9(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent9 = params.trigEvent9/sum(trials.reward);
        
            
            % make matrix
            RL_event = concat_event(params);
            
            % when align pupil signal to cue
            params.trigTime = trialData.cueTimes;
            params.xtitle = 'Time from cue (s)';
            params.window = [-3:0.1:5];
            params.nback = 0;       %how many trials back to regress against
            params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            params.interaction = false;
            params.ifplot = 0;
            %only perform analysis on this subset of trials
            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
            tlabel={'C(n)','C(n-1)','R(n-1)','deltaQ','RPE', 'deltaK','CKE','Reward Rate', 'Cumulavtive reward'};
            
            %reg_cr=linear_regr( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params );
            %         reg_cr_ctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
            %         save(saveMLRmatpath, 'reg_cr','reg_cr_ctrl');
            
            params.xtitle = {'Time from cue (s)'};
            LR_t = cells.t;
            LR_dFF = cells.dFF;
            parfor j=1:numel(LR_dFF )
                if length(LR_t) > length(LR_dFF {1})
                    reg_cr{j}=linear_regr( LR_dFF {j}, LR_t(1:length(LR_dFF {1})), RL_event, params.trigTime, trialMask, params );
                else
                    reg_cr{j}=linear_regr( LR_dFF {j}, LR_t, RL_event, params.trigTime, trialMask, params );
                end
            end
            MP_plot_regr(reg_cr,[],params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','MLR-valueupdating');    %png format
            saveas(gcf, 'MLR-valueupdating', 'fig');
            
            %% running control multilinear regression
            % shuffle every factor one by one, keeping other factors intact
            
            % check
            
            %             % construct every regression factor
            %             params_ctrl = params;
            %
            %             % concatenate it into a matrix
            %             params_ctrlMat = [];
            %             fields = fieldnames(params_ctrl);
            %             for jj = 1:length(fields)
            %                 if contains(fields{jj},'trigEvent')
            %                     params_ctrlMat = [params_ctrlMat params_ctrl.(fields{jj})];
            %                 end
            %             end
            %
            %             % shul
            %             fieldname={'go'};
            %             trialMask = getMask(trials,fieldname);
            %
            %
            %             params.xtitle = 'Time from cue (s)';
            %             params.window = [-3:0.1:5];
            %             params.nback = 0;       %how many trials back to regress against
            %             params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            %             params.interaction = false;
            %             %params.trigTime = trialData.cueTimes;
            %             params.ifplot = 0;
            %             %only perform analysis on this subset of trials
            %
            %             % iterate through all 9 factors, shuffle name one by one to get
            %             % the control regression for every factor
            %             tlabel={'C(n)','R(n)','C(n)xR(n)','C(n-1)','R(n-1)', 'C(n-1)xR(n-1)','QL-QR', 'ChosenQ', 'Reward Rate', 'Cumulavtive reward'};
            %
            %
            %             reg_cr_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, params_ctrlMat, params.trigTime, trialMask, params, tlabel);
            %             reg_cr_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, params_ctrlMat,  params.trigTime, trialMask, params, tlabel );
            %
            %             save(saveMLRmatpath, 'reg_cr','reg_cr_ctrl');
            save(saveRegName, 'reg_cr');
            close all;
        else
            display('Regression already done');
        end
        
    end
end
end