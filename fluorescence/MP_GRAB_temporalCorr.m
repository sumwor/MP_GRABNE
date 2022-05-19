function MP_GRAB_temporalCorr(dataIndex)

% calculate temporal correlation of different varible coefficients within
% grids

% check if coefficient sign changes over time
% for eachgrid, calculate the number of + time point and - time point
% then average

nFiles = size(dataIndex,1);

for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    %load(fullfile(fn_beh.folder,fn_beh.name));


    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1

        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end

        cd(savefluofigpath);


        %load(fullfile(fn_fluo.folder,fn_fluo.name));

        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end



           %% set threshold
        Thresh.pos = 0.2; Thresh.neg = -0.2; % threshold to find area with pos/neg modulation
        sigThresh.value = 0.05; sigThresh.alpha = 0.01;  % threshold to find area with significant modulation

        tic
    
        %% check other variables using regression results (should focus on significance rather than coefficient)

        % variable list:
        % regression 1: cn,rn,cn+1,cn-1,rn-1,cn*rn,average_r,cumulative_r
        % regression 2: dQ,dK,chosenQ,chosenK
        % regression 3: RPE, CKE,
        saveregmaskpath = fullfile(savematpath,'regressionTempCorr.mat');  % file path to save the results

        if exist(saveregmaskpath)
            load(saveregmaskpath);
        end

        % need to get the boundaries of possible coefficient value of all
        % regressions
        saveCRName = fullfile(savematpath,'regCR_norm.mat');
         if exist(saveCRName)
            reg_cr{1} = load(saveCRName);
         end
        saveRLName = fullfile(savematpath,'regRL_norm.mat');
          if exist(saveRLName)
            reg_cr{2} = load(saveRLName);
         end
        saveRPEName = fullfile(savematpath,'regRL_RPE_norm.mat');
          if exist(saveRPEName)
            reg_cr{3} = load(saveRPEName);
          end
        
         % get corresponding coefficient
         numGrids = numel(reg_cr{1}.reg_cr);
          varName = {'cn+1','cn','cn-1','rn','rn-1','xn','xn-1','ave_r','Cum_r','dQ','chosenQ','dK','chosenK','RPE','CKE'};
        varInd = [2,3,4,7,8,11,12,14,15,8,9,10,11,6,8];
        varReg = [1,1,1,1,1,1,1,1,1,2,2,2,2,3,3];
        zerolag = length(reg_cr{1}.reg_cr{1}.regr_time);
       
       

        tempCorrMat = cell(length(varInd),length(varInd));
        corrCoeff0 = zeros(length(varInd),length(varInd),numGrids);
        corrMax = zeros(length(varInd),length(varInd),numGrids);
        corrLag = zeros(length(varInd),length(varInd),numGrids);
        % calculate crosscorrelation
        %if ~exist(saveregVarcorrpath)
        for uu = 1:length(varName)
            for vv=uu:length(varName)

                tempCorrMat{uu,vv} = getRegTempCorr(reg_cr{varReg(uu)},reg_cr{varReg(vv)},varName{uu},varName{vv},varInd(uu),varInd(vv), sigThresh, savefluofigpath);
                corrCoeff0(uu,vv,:) = squeeze(tempCorrMat{uu,vv}(zerolag,1,:));
                [corrMax(uu,vv,:),maxInd] = max(squeeze(tempCorrMat{uu,vv}(:,1,:)),[],1);
                lagVec = squeeze(tempCorrMat{uu,vv}(:,2,:));
                corrLag(uu,vv,:) = lagVec(maxInd);
            end
        %end
        end
        
        % save the data
        save(saveregmaskpath,'tempCorrMat','corrCoeff0','corrMax','corrLag');

%% with in variables, calculate temporal correlation with a same reference grid
    saveregpath = fullfile(savematpath,'regressionVarTempCorr.mat');  % file path to save the results

        %% regression 1 -------------------------------------------------------------------------
             % get the choice from linear regression
        %if ~exist('choicetempData','var') % if choice regression mask not computed
            label = 'choice';
            choiceInd = 3;

            choicetempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,choiceInd,sigThresh,savefluofigpath);
        %end

        % outcome regression mask
        
        %if ~exist('outcometempData','var') % if choice regression mask not computed
            label = 'outcome';
            outcomeInd = 7;
            outcometempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,outcomeInd,sigThresh,savefluofigpath);
        %end
        % save the result

         % c(n+1) regression mask
        %if ~exist('cn_1tempData','var') % if choice regression mask not computed
            label = 'upcoming choice';
            cn_1Ind = 2;
            cn_1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,cn_1Ind,sigThresh,savefluofigpath);
        %end

        % r(n+1) regression mask
        %if ~exist('rn_1tempData','var') % if choice regression mask not computed
            label = 'upcoming reward';
            rn_1Ind = 6;
            rn_1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,rn_1Ind,sigThresh,savefluofigpath);
        %end

        % c(n-1) regression mask
        %if ~exist('cn__1tempData','var') % if choice regression mask not computed
            label = 'previous choice';
            cn__1Ind = 4;
            cn__1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,cn__1Ind,sigThresh,savefluofigpath);
        %end

        % r(n-1) regression mask
        %if ~exist('rn__1tempData','var') % if choice regression mask not computed
            label = 'previous reward';
            rn__1Ind = 8;
            rn__1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,rn__1Ind,sigThresh,savefluofigpath);
        %end

        % x(n) regression mask
        %if ~exist('xntempData','var') % if choice regression mask not computed
            label = 'interaction';
            xnInd = 11;
            xntempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,xnInd,sigThresh,savefluofigpath);
        %end

         %if ~exist('xn__1tempData','var')
         label = 'previous interaction';
            xnInd = 12;
            xn__1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,xnInd,sigThresh,savefluofigpath);
         %end

        % average reward regression mask
        %if ~exist('ave_rtempData','var') % if choice regression mask not computed
            label = 'average reward';
            ave_rInd = 14;
            ave_rtempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,ave_rInd,sigThresh,savefluofigpath) ;
        %end

        % cumulative reward regression mask
        %if ~exist('cum_rtempData','var') % if choice regression mask not computed
            label = 'cumulative reward';
            cum_rInd = 15;
            cum_rtempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,cum_rInd,sigThresh,savefluofigpath) ;
        %end

        %% regression 2-----------------------------------------------------------------------
        %if ~exist('dQtempData','var') % if choice regression mask not computed
            label = 'delta Q';
            dQInd = 8;
            dQtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,dQInd,sigThresh,savefluofigpath) ;
        %end

        %if ~exist('chosenQtempData','var') % if choice regression mask not computed
            label = 'chosen Q';
            chosenQInd = 9;
            chosenQtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,chosenQInd, sigThresh,savefluofigpath) ;
        %end

        %if ~exist('dKtempData','var') % if choice regression mask not computed
            label = 'delta K';
            dKInd = 10;
            dKtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,dKInd,sigThresh,savefluofigpath) ;
        %end

        %if ~exist('chosenKtempData','var') % if choice regression mask not computed
            label = 'chosen K';
            chosenKInd = 11;
            chosenKtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,chosenKInd, sigThresh,savefluofigpath) ;
        %end

        %% regression 3-----------------------------------------------------------------------

        %if ~exist('RPEtempData','var') % if choice regression mask not computed
            label = 'RPE';
            RPEInd = 6;
            RPEtempData = getRegautoCorrData(reg_cr{3}.reg_cr,label,RPEInd,sigThresh,savefluofigpath) ;
        %end

        %if ~exist('CKEtempData','var') % if choice regression mask not computed
            label = 'CKE';
            CKEInd = 8;
            CKEtempData = getRegautoCorrData(reg_cr{3}.reg_cr,label,CKEInd,sigThresh,savefluofigpath) ;
        %end

        %% save the results
        save(saveregpath,'outcometempData','choicetempData','cn_1tempData','cn__1tempData',...
                             'rn_1tempData','rn__1tempData', 'xntempData', 'xn__1tempData', 'ave_rtempData', 'cum_rtempData',...
                             'dQtempData', 'chosenQtempData','dKtempData','chosenKtempData',...
                             'RPEtempData', 'CKEtempData')
        
       
    end

        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
    toc
end
end


