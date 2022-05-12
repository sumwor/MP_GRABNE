function MP_GRAB_clustering(dataIndex);

% try: hierachical clustering; spectral clustering


nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
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
        
        warning('off','all')
        
           %% set threshold
        Thresh.pos = 0.2; Thresh.neg = -0.2; % threshold to find area with pos/neg modulation
        sigThresh.value = 0.05; sigThresh.alpha = 0.01;  % threshold to find area with significant modulation
        
        tic
        
        %% check other variables using regression results (should focus on significance rather than coefficient)
        
        % variable list:
        % regression 1: cn,rn,cn+1,cn-1,rn-1,cn*rn,average_r,cumulative_r
        % regression 2: dQ,dK,chosenQ,chosenK
        % regression 3: RPE, CKE,
        saveregclusterpath = fullfile(savematpath,'cluster.mat');  % file path to save the results
        
        if exist(saveregclusterpath)
            load(saveregclusterpath);
        end
        
        % need to get the boundaries of possible coefficient value of all
        % regressions
        saveCRName = fullfile(savematpath,'regCR_norm.mat');
         if exist(saveCRName)
            reg1 = load(saveCRName);
         end
        saveRLName = fullfile(savematpath,'regRL_norm.mat');
          if exist(saveRLName)
            reg2 = load(saveRLName);
         end
        saveRPEName = fullfile(savematpath,'regRL_RPE_norm.mat');
          if exist(saveRPEName)
            reg3 = load(saveRPEName);
          end
         
          % data should be save in 28*28*n*80
          % or 748*n*80?
        clusterData = zeros(28,28,9,80); 
        % using choice,reward,interaction,aveR,cumR,dQ,dK,RPE,CKE for now
        %% regression 1 -------------------------------------------------------------------------
             % get the choice from linear regression
        %if ~exist('choiceRegData','var') % if choice regression mask not computed
            label = 'choice';
            choiceInd = 3;
            choice = selectRg2D(reg1.reg_cr,choiceInd,[]);
            sigChoice=choice.coeff;
            sigChoice(choice.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,1,:) = sigChoice;
         
        %end
        
        % outcome regression mask
       % if ~exist('outcomeRegData','var') % if choice regression mask not computed
            label = 'outcome';
            outcomeInd = 7;
            outcome = selectRg2D(reg1.reg_cr,outcomeInd,[]);
            sigOutcome=outcome.coeff;
            sigOutcome(outcome.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,2,:) = sigOutcome;
        
        % x(n) regression mask
        %if ~exist('xnRegData','var') % if choice regression mask not computed
            label = 'interaction';
            xnInd = 11;
             interaction = selectRg2D(reg1.reg_cr,xnInd,[]);
            siginteraction=interaction.coeff;
            siginteraction(interaction.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,3,:) = siginteraction;
        %end
        
        % average reward regression mask
        %if ~exist('ave_rRegData','var') % if choice regression mask not computed
            label = 'average reward';
            ave_rInd = 14;
            ave_r = selectRg2D(reg1.reg_cr,ave_rInd,[]);
            sigave_r=ave_r.coeff;
            sigave_r(ave_r.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,4,:) = sigave_r;
        %end
        
        % cumulative reward regression mask
        %if ~exist('cum_rRegData','var') % if choice regression mask not computed
            label = 'cumulative reward';
            cum_rInd = 15;
            cum_r = selectRg2D(reg1.reg_cr,cum_rInd,[]);
            sigcum_r=cum_r.coeff;
            sigcum_r(cum_r.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,5,:) = sigcum_r;
        %end
        
        %% regression 2-----------------------------------------------------------------------
        
        %if ~exist('dQRegData','var') % if choice regression mask not computed
            label = 'delta Q';
            dQInd = 8;
            dQ = selectRg2D(reg2.reg_cr,dQInd,[]);
            sigdQ=dQ.coeff;
            sigdQ(dQ.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,6,:) = sigdQ;
    %end
        %if ~exist('chosenQRegData','var') % if choice regression mask not computed

        %if ~exist('dKRegData','var') % if choice regression mask not computed
            label = 'delta K';
            dKInd = 10;
            dK = selectRg2D(reg2.reg_cr,dKInd,[]);
            sigdK=dK.coeff;
            sigdK(dK.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,7,:) = sigdK;
        %end
        
        
        %% regression 3-----------------------------------------------------------------------
        
        %if ~exist('RPERegData','var') % if choice regression mask not computed
            label = 'RPE';
            RPEInd = 6;
            RPE = selectRg2D(reg3.reg_cr,RPEInd,[]);
            sigRPE=RPE.coeff;
            sigRPE(RPE.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,8,:) = sigRPE;
        %end
        
        %if ~exist('CKERegData','var') % if choice regression mask not computed
            label = 'CKE';
            CKEInd = 8;
            CKE = selectRg2D(reg3.reg_cr,CKEInd,[]);
            sigCKE=CKE.coeff;
            sigRPE(CKE.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,9,:) = sigCKE;
        %end
        
        HCluster = reshape(clusterData,784,720);
        T = clusterdata(HCluster,'Linkage','ward','SaveMemory','on','Maxclust',20);
        
        % plot average coefficient for different clusters
        ClusterMat = zeros(28,28);
        for tt = 1:length(T)
            ind1 = floor((tt-1)/28)+1;ind2 =mod(tt,28);
            if ind2==0
                ind2 = 28;
            end
            ClusterMat(ind1,ind2) = T(tt);
        end
        figure;imagesc(ClusterMat)
        
        tree = linkage(HCluster,'average');

        figure;dendrogram(tree);
        
        rng('default');  % For reproducibility
X = rand(20000,3);
 T = clusterdata(X,'Linkage','ward','SaveMemory','on','Maxclust',4);
        
        %% save the results
        save(saveregmaskpath,'outcomeRegData','choiceRegData','cn_1RegData','cn__1RegData',...
                             'rn__1RegData', 'xnRegData', 'ave_rRegData', 'cum_rRegData',...
                             'dQRegData', 'chosenQRegData','dKRegData','chosenKRegData',...
                             'RPERegData', 'CKERegData')
    end
    
        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
    toc
end
end