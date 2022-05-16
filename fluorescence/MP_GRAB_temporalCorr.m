function MP_GRAB_temporalCorr(dataIndex)

% calculate temporal correlation of different varible coefficients within
% grids

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
        
        % correlation clustering
        
%         colors=cbrewer('div','RdBu',256);
%         colors=flipud(colors);
%         colorRange = [-1 1];
% 
%         figure;image(nanmean(corrMax,3),'CDataMapping','scaled')
%          axis square;
%     colormap(colors);
%     caxis([colorRange(1) colorRange(2)]);
%         
        
          save(saveregmaskpath,'tempCorrMat','corrCoeff0','corrMax','corrLag');
    end

        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
    toc
end
end


