function MP_GRAB_selectivitySummary(dataIndex, save_path_fluo)

%% summarize selectivity results, as well as cross correlation results
nFiles = size(dataIndex,1);

%% selectivity area
% RPE/R
RPE_R_sig = zeros(1, nFiles);RPE_sig = zeros(1, nFiles);R_sig = zeros(1,nFiles);
sigVar = zeros(28,28,80,nFiles,15);  % load significant grids for 15 variables

% cn/xn
% cn_xn_sig = zeros(1,nFiles); ncn_xn_sig = zeros(1,nFiles);
% rn_xn_sig = zeros(1,nFiles); nrn_xn_sig = zeros(1,nFiles);
% cnrn_xn_sig = zeros(1,nFiles); ncnrn_xn_sig = zeros(1,nFiles);
cnSigList = [];
rnSigList = [];
xnSigList = [];
posRPESigList = [];
negRPESigList = [];
% get the pos/neg coefficient percentage data
choicePN = [];outcomePN = [];RPEPN=[];xnPN=[];cn__1PN=[];rn__1PN=[];xn__1PN = [];
dQPN = [];chosenQPN = []; dKPN = [];chosenKPN = [];CKEPN = [];


for ii = 1:nFiles
    tic
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));


    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1

        % make folders to save analysis and plots
        warning('off','all')
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        %% linear regression selection
        saveregmaskpath = fullfile(savematpath,'regressionMask.mat');
        fn_regSel = dir(saveregmaskpath);
        if length(fn_regSel) == 1
            load(saveregmaskpath);
        else
            display('Selectivity not computed');
        end


        % RPE/reward compare
        % check how many grids are both signficant/RPE sig only/reward sig
        % only
        % get sum of Mask first
%         sigVar(:,:,:,ii,1) = cn_1RegData.sigMask;
%         sigVar(:,:,:,ii,2) = choiceRegData.sigMask;
%         sigVar(:,:,:,ii,3) = cn__1RegData.sigMask;
%         sigVar(:,:,:,ii,4) = outcomeRegData.sigMask;
%        sigVar(:,:,:,ii,5) = rn__1RegData.sigMask;
%        sigVar(:,:,:,ii,6) = xnRegData.sigMask;
%        sigVar(:,:,:,ii,7) = xn__1RegData.sigMask;
%        sigVar(:,:,:,ii,8) = ave_rRegData.sigMask;
%        sigVar(:,:,:,ii,9) = cum_rRegData.sigMask;
%        sigVar(:,:,:,ii,10) = dQRegData.sigMask;
%           sigVar(:,:,:,ii,11) = chosenQRegData.sigMask;
%              sigVar(:,:,:,ii,12) = dKRegData.sigMask;
%                 sigVar(:,:,:,ii,13) = chosenKRegData.sigMask;
%        sigVar(:,:,:,ii,14) = RPERegData.sigMask;
%        sigVar(:,:,:,ii,15) = CKERegData.sigMask;

%         RPE_R_sig(ii) = sum(sum((~RPEnotSig) & (~RnotSig)));
%         RPE_sig(ii) = sum(sum((~RPEnotSig)&(RnotSig)));
%         R_sig(ii) = sum(sum((RPEnotSig)&(~RnotSig)));
        % get significant numbers of interaction and choice
        % cn/rn/xn
%         sumGrid = sum(choiceRegData.sigGrid(:)|xnRegData.sigGrid(:)|outcomeRegData.sigGrid(:));
%         cn_xn_sig(ii) = sum(choiceRegData.sigGrid(:)&xnRegData.sigGrid(:))/sumGrid;
%         ncn_xn_sig(ii) = sum(~choiceRegData.sigGrid(:)&xnRegData.sigGrid(:))/sumGrid;
%         rn_xn_sig(ii) = sum(xnRegData.sigGrid(:)&outcomeRegData.sigGrid(:))/sumGrid;
%         nrn_xn_sig(ii) = sum(xnRegData.sigGrid(:)&~outcomeRegData.sigGrid(:))/sumGrid;
%         cnrn_xn_sig(ii) = sum(choiceRegData.sigGrid(:)&xnRegData.sigGrid(:)&outcomeRegData.sigGrid(:))/sumGrid;
%         ncnrn_xn_sig(ii) = sum(~choiceRegData.sigGrid(:)&xnRegData.sigGrid(:)&~outcomeRegData.sigGrid(:))/sumGrid;
        cnSigList = [cnSigList,choiceRegData.sigGrid(:)'];
        rnSigList = [rnSigList,outcomeRegData.sigGrid(:)'];
        xnSigList = [xnSigList,xnRegData.sigGrid(:)'];
        posRPESigList = [posRPESigList, posRPERegData.sigGrid(:)'];
        negRPESigList = [negRPESigList, negRPERegData.sigGrid(:)'];

        % get pos/neg coefficient
%         choicePN = cat(4,choicePN,choiceRegData.coeffSignSig);
%          outcomePN = cat(4,outcomePN,outcomeRegData.coeffSignSig);
%          RPEPN = cat(4,RPEPN,RPERegData.coeffSignSig);
%          xnPN=cat(4,xnPN,xnRegData.coeffSignSig);
%          cn__1PN=cat(4,cn__1PN,cn__1RegData.coeffSignSig);
%          rn__1PN=cat(4,rn__1PN,rn__1RegData.coeffSignSig);
%          xn__1PN = cat(4,xn__1PN,xn__1RegData.coeffSignSig);
%          dQPN = cat(4,dQPN,dQRegData.coeffSignSig);
%          chosenQPN = cat(4,chosenQPN,chosenQRegData.coeffSignSig);
%          dKPN = cat(4,dKPN,dKRegData.coeffSignSig);
%          chosenKPN = cat(4,chosenKPN,chosenKRegData.coeffSignSig);
%          CKEPN = cat(4,CKEPN,CKERegData.coeffSignSig);

    end
    toc
end

% reshape
% s1 = size(choicePN,1)*size(choicePN,2);s2 = size(choicePN,3);s3 = size(choicePN,4);
% choicePN = reshape(choicePN,s1,s2,s3);
% outcomePN = reshape(outcomePN,s1,s2,s3);
% RPEPN=reshape(RPEPN,s1,s2,s3);
% xnPN=reshape(xnPN,s1,s2,s3);
% cn__1PN=reshape(cn__1PN,s1,s2,s3);
% rn__1PN=reshape(rn__1PN,s1,s2,s3);
% xn__1PN = reshape(xn__1PN,s1,s2,s3);
% dQPN = reshape(dQPN,s1,s2,s3);
% chosenQPN =reshape(chosenQPN,s1,s2,s3);
% dKPN = reshape(dKPN,s1,s2,s3);
% chosenKPN = reshape(chosenKPN,s1,s2,s3);
% CKEPN = reshape(CKEPN,s1,s2,s3);

% plot overall average coefficient
% choice/outcome/interaction/RPE
% s1 = size(choicePN,1)*size(choicePN,3);
% posVar1 = [reshape(choicePN(:,1,:),1,s1);reshape(outcomePN(:,1,:),1,s1);reshape(xnPN(:,1,:),1,s1);reshape(RPEPN(:,1,:),1,s1)];
% negVar1 = [reshape(choicePN(:,2,:),1,s1);reshape(outcomePN(:,2,:),1,s1);reshape(xnPN(:,2,:),1,s1);reshape(RPEPN(:,2,:),1,s1)];
% violin plot
%figure
%v=violinplot(posVar1',animalIndex,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0],'ShowData',false);
%  colors = [235, 83, 83;249, 217, 35;54, 174, 124;24, 116, 152]/255;
% figure;violinplot(posVar1',{'choice','outcome','interaction','RPE'},'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylabel('Fraction of positive time');
% print(gcf,'-dpng',fullfile(save_path_fluo,'Fraction of positive time -choice'));
% saveas(gcf, fullfile(save_path_fluo,'Fraction of positive time -choice'), 'fig');
% %
% posVar2 = [reshape(dQPN(:,1,:),1,s1);reshape(chosenQPN(:,1,:),1,s1);reshape(dKPN(:,1,:),1,s1);reshape(chosenKPN(:,1,:),1,s1)];
% negVar2 = [reshape(dQPN(:,2,:),1,s1);reshape(chosenQPN(:,2,:),1,s1);reshape(dKPN(:,2,:),1,s1);reshape(chosenKPN(:,2,:),1,s1)];
% figure;violinplot(posVar2',{'dQ','chosenQ','dK','chosenK'},'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylabel('Fraction of positive time');
% print(gcf,'-dpng',fullfile(save_path_fluo,'Fraction of positive time -latent'));
% saveas(gcf, fullfile(save_path_fluo,'Fraction of positive time -latent'), 'fig');

% individual animals
% animalList = unique(dataIndex.Animal);
% for aa = 1:numel(animalList)
%     sessionInclude = [];
%     for tt = 1:size(dataIndex,1)
%         if dataIndex.Animal{tt}==animalList{aa}
%             sessionInclude = [sessionInclude, tt];
%         end
%     end
%     s1 = size(choicePN,1)*length(sessionInclude);
%     posVar1 = [reshape(choicePN(:,1,sessionInclude),1,s1);reshape(outcomePN(:,1,sessionInclude),1,s1);reshape(xnPN(:,1,sessionInclude),1,s1);reshape(RPEPN(:,1,sessionInclude),1,s1)];
%     negVar1 = [reshape(choicePN(:,2,sessionInclude),1,s1);reshape(outcomePN(:,2,sessionInclude),1,s1);reshape(xnPN(:,2,sessionInclude),1,s1);reshape(RPEPN(:,2,sessionInclude),1,s1)];
%
%     figure;violinplot(posVar1',{'choice','outcome','interaction','RPE'},'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylabel('Fraction of positive time');
% title(animalList{aa});
% print(gcf,'-dpng',fullfile(save_path_fluo,['Fraction of positive time -choice ', animalList{aa}]));
% saveas(gcf, fullfile(save_path_fluo,['Fraction of positive time -choice ', animalList{aa}]), 'fig');

% posVar2 = [reshape(dQPN(:,1,sessionInclude),1,s1);reshape(chosenQPN(:,1,sessionInclude),1,s1);reshape(dKPN(:,1,sessionInclude),1,s1);reshape(chosenKPN(:,1,sessionInclude),1,s1)];
% negVar2 = [reshape(dQPN(:,2,sessionInclude),1,s1);reshape(chosenQPN(:,2,sessionInclude),1,s1);reshape(dKPN(:,2,sessionInclude),1,s1);reshape(chosenKPN(:,2,sessionInclude),1,s1)];
% figure;violinplot(posVar2',{'dQ','chosenQ','dK','chosenK'},'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylabel('Fraction of positive time');
% title(animalList{aa});
% print(gcf,'-dpng',fullfile(save_path_fluo,['Fraction of positive time -latent', animalList{aa}]));
% saveas(gcf, fullfile(save_path_fluo,['Fraction of positive time -latent', animalList{aa}]), 'fig');

% end

% close all;

% figure;violinplot(negVar1',{'choice','outcome','interaction','RPE'});
% figure;violinplot(posVar1'./negVar1',{'choice','outcome','interaction','RPE'});
% bar plot

% check the interactionn with choice/outcome

RInd = 1:length(cnSigList);
cnInd = RInd(cnSigList==1);
rnInd = RInd(rnSigList==1);
xnInd = RInd(xnSigList==1);
posRPEInd = RInd(posRPESigList==1);
negRPEInd = RInd(negRPESigList==1);


label = {'choice','outcome','interaction'};
vennPlot(cnInd,rnInd,xnInd,label,save_path_fluo);

label = {'Outcome','posRPE','negRPE'};
vennPlot(rnInd,posRPEInd,negRPEInd,label,save_path_fluo);

h = vennEulerDiagram({cnInd;rnInd;xnInd}, 'drawProportional', true, 'SetLabels', ["choice"; "reward"; "interaction"]);
title('Ach c/r/x significant ROIs');
print(gcf,'-dpng',fullfile(save_path_fluo,'Ratio of significant grids (c-r-x)-venn'));
saveas(gcf, fullfile(save_path_fluo,'Ratio of significant grids (c-r-x)-venn'), 'fig');

data = [cn_xn_sig;ncn_xn_sig;rn_xn_sig;nrn_xn_sig;cnrn_xn_sig;ncnrn_xn_sig];
figure
violinplot(data',{'c&x','nc&x','r&x','nr&x','cr&x','ncr&x'});
ylim([0 1]);
title('Ratio of significant grids');
print(gcf,'-dpng',fullfile(save_path_fluo,'Ratio of significant grids (c-r-x)'));
saveas(gcf, fullfile(save_path_fluo,'Ratio of significant grids (c-r-x)'), 'fig');


% cats={'RPE & R','RPE only','R only'};
% figure;
% violinplot([RPE_R_sig;RPE_sig;R_sig]',cats);
% ylabel('Number of significant grids');
%
% print(gcf,'-dpng',fullfile(save_path_fluo,'Sig_RPE_R'));    %png format
% saveas(gcf, fullfile(save_path_fluo,'Sig_RPE_R'), 'fig');
% % paired t test?
% [h1,p1] = ttest(RPE_sig,R_sig)

%% cross correlation
% predictors of interest: RPE/Reward/Cumulative reward/average reward
% choice
% deltaQ/deltaK

choice_sig_dis = [];choice_center_coeff = [];choice_regPval = [];%choice_sigmax_ind = [];choice_sigmax_val=[];
outcome_sig_dis = [];outcome_center_coeff = [];outcome_regPval = [];%outcome_sigmax_ind = [];outcome_sigmax_val=[];
rn_1_sig_dis  = [];rn_1_center_coeff=[];rn_1_regPval = [];
ave_r_sig_dis = [];ave_r_center_coeff = [];ave_r_regPval = [];
cum_r_sig_dis = [];cum_r_center_coeff = [];cum_r_regPval = [];
dQ_sig_dis = [];dQ_center_coeff = [];dQ_regPval = [];
dK_sig_dis = [];dK_center_coeff = [];dK_regPval = []; dK_sigCorr = [];
chosenQ_sig_dis = [];chosenQ_center_coeff = [];chosenQ_regPval = [];
chosenK_sig_dis = [];chosenK_center_coeff = [];chosenK_regPval = [];
RPE_sig_dis = [];RPE_center_coeff = [];RPE_regPval = [];%RPE_sigmax_ind = [];RPE_sigmax_val=[];


%% crosscorrelation
xCorrMat = zeros(15,15,80,nFiles);
xSigMat = zeros(15,15,80,nFiles);
regCoeffMat = zeros(15,80,nFiles); % save the regression coefficient

% regression results


for ii = 1:nFiles
    tic
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));


    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1

        % make folders to save analysis and plots
        warning('off','all')
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        %% linear regression selection
        saveregXcorrpath = fullfile(savematpath,'regressionXcorr.mat');
        if exist(saveregXcorrpath)
            load(saveregXcorrpath);
        end
        % choice
        choice_center_coeff=[choice_center_coeff;squeeze(choiceRegData.sigCorr(28,28,:))'];
        choice_sig_dis = [choice_sig_dis;choiceRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        choice_regPval = cat(3,choice_regPval,choiceRegData.regPval);
        %
         % outcome
        outcome_center_coeff=[outcome_center_coeff;squeeze(outcomeRegData.sigCorr(28,28,:))'];
        outcome_sig_dis = [outcome_sig_dis;outcomeRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        outcome_regPval = cat(3,outcome_regPval,outcomeRegData.regPval);

         % rn+1
        rn_1_center_coeff=[rn_1_center_coeff;squeeze(rn_1RegData.sigCorr(28,28,:))'];
        rn_1_sig_dis = [rn_1_sig_dis;rn_1RegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        rn_1_regPval = cat(3,rn_1_regPval,rn_1RegData.regPval);

        % average reward
        ave_r_center_coeff=[ave_r_center_coeff;squeeze(ave_rRegData.sigCorr(28,28,:))'];
        ave_r_sig_dis = [ave_r_sig_dis;ave_rRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        ave_r_regPval = cat(3,ave_r_regPval,ave_rRegData.regPval);

        % cumulative reward
        cum_r_center_coeff=[cum_r_center_coeff;squeeze(cum_rRegData.sigCorr(28,28,:))'];
        cum_r_sig_dis = [cum_r_sig_dis;cum_rRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        cum_r_regPval = cat(3,cum_r_regPval,cum_rRegData.regPval);

        % dQ
        dQ_center_coeff=[dQ_center_coeff;squeeze(dQRegData.sigCorr(28,28,:))'];
        dQ_sig_dis = [dQ_sig_dis;dQRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        dQ_regPval = cat(3,dQ_regPval,dQRegData.regPval);

        % dK
        dK_center_coeff=[dK_center_coeff;squeeze(dKRegData.sigCorr(28,28,:))'];
        dK_sig_dis = [dK_sig_dis;dKRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        dK_regPval = cat(3,dK_regPval,dKRegData.regPval);
        dK_sigCorr = cat(4,dK_sigCorr,dKRegData.sigCorr);

        % chosenQ
        chosenQ_center_coeff=[chosenQ_center_coeff;squeeze(chosenQRegData.sigCorr(28,28,:))'];
        chosenQ_sig_dis = [chosenQ_sig_dis;chosenQRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        chosenQ_regPval = cat(3,chosenQ_regPval,chosenQRegData.regPval);

        % chosenK
        chosenK_center_coeff=[chosenK_center_coeff;squeeze(chosenKRegData.sigCorr(28,28,:))'];
        chosenK_sig_dis = [chosenK_sig_dis;chosenKRegData.sigMaxDis];
        %choice_sigmax_ind =  cat(3,choice_sigmax_ind,choiceRegData.sigMaxInd);
        %choice_sigmax_val = [choice_sigmax_val;choiceRegData.sigMaxValue];
        chosenK_regPval = cat(3,chosenK_regPval,chosenKRegData.regPval);

        % RPE
        RPE_center_coeff=[RPE_center_coeff;squeeze(RPERegData.sigCorr(28,28,:))'];
        RPE_sig_dis = [RPE_sig_dis;RPERegData.sigMaxDis];
%         RPE_sigmax_ind =  cat(3,RPE_sigmax_ind,RPERegData.sigMaxInd);
%         RPE_sigmax_val = [RPE_sigmax_val;RPERegData.sigMaxValue];
        RPE_regPval = cat(3,RPE_regPval,RPERegData.regPval);

        % pos/neg coefficient ratio
        %% cross correlation results
         saveregVarcorrpath = fullfile(savematpath,'regressionVarcorr.mat');
         load(saveregVarcorrpath);
         for uu = 1:size(xCorrCell,1)
             for vv = uu:size(xCorrCell,2)
                 xCorrMat(uu,vv,:,ii) = xCorrCell{uu,vv}.coeffCorr(28,28,:);
                 xSigMat(uu,vv,:,ii) = xCorrCell{uu,vv}.sigCorr(28,28,:);
             end
         end
    end
end

%% check the cross-correlation plot

%% ---------------choice
label = 'choice significance';
maxValue.coeff= choice_center_coeff;
offset.coeff = choice_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = choice_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);
% check the coordinate of maximum correlation

%% ---------------outcome
label = 'outcome significance';
maxValue.coeff= outcome_center_coeff;
offset.coeff = outcome_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = outcome_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);

%% ---------------upcoming reward
label = 'rn+1 significance';
maxValue.coeff= rn_1_center_coeff;
offset.coeff = rn_1_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = rn_1_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);

%% ---------------average reward
label = 'average reward significance';
maxValue.coeff= ave_r_center_coeff;
offset.coeff = ave_r_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = ave_r_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);
% check the coordinate of maximum coefficient between 0-2s

%% ---------------cumulative reward
label = 'cumulative reward significance';
maxValue.coeff= cum_r_center_coeff;
offset.coeff = cum_r_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = cum_r_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);
% check t=0.05

%% ---------------dQ
label = 'dQ significance';
maxValue.coeff= dQ_center_coeff;
offset.coeff = dQ_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = dQ_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);

%% ---------------dK
label = 'dK significance';
maxValue.coeff= dK_center_coeff;
offset.coeff = dK_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = dK_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);
% average crosscorrelation?

%% ---------------chosen Q
label = 'chosenQ significance';
maxValue.coeff= chosenQ_center_coeff;
offset.coeff = chosenQ_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = chosenQ_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);

%% ---------------chosen K
label = 'chosenK significance';
maxValue.coeff= chosenK_center_coeff;
offset.coeff = chosenK_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = chosenK_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);

%% ---------------RPE
% plot average
label = 'RPE significance';
maxValue.coeff= RPE_center_coeff;
offset.coeff = RPE_sig_dis;
maxValue.t = -2.95:0.1:4.85; offset.t = -2.95:0.1:4.85;
pValue.coeff = RPE_regPval;pValue.t = -2.95:0.1:4.95;
thresh = 0.01;
plot_xcorrCoeff(pValue,maxValue,offset,thresh,save_path_fluo,label);



% plot the coordinate where xcorr reaches maximum
% figure;
%
% scatter(RPE_sigmax_ind(1,:,1),RPE_sigmax_ind(2,:,1),40,1:size(RPE_sigmax_ind,2),'filled');
% colormap('jet')
% colorbar

%% summarize the crosscorrelation results
% bootstrap the results

xCorrBoot = cell(size(xCorrMat,1),size(xCorrMat,1));
xCorrBootAve =  zeros(size(xCorrMat,1),size(xCorrMat,1),size(xCorrMat,3));
xCorrBootPval = ones(size(xCorrMat,1),size(xCorrMat,1),size(xCorrMat,3));


xSigBoot = cell(size(xCorrMat,1),size(xCorrMat,1));
xSigBootAve =  zeros(size(xCorrMat,1),size(xCorrMat,1),size(xCorrMat,3));
xSigBootPval = ones(size(xCorrMat,1),size(xCorrMat,1),size(xCorrMat,3));

for tt = 1:size(xCorrMat,1)
    for rr = tt:size(xCorrMat,2)
        if tt~= rr
            bootData.coeff = squeeze(xCorrMat(tt,rr,:,:))';
            xCorrBoot{tt,rr} = getBootstrp(bootData,0,0.05);
            xCorrBootAve(tt,rr,:) = xCorrBoot{tt,rr}.coeff_bootave;

            bootSigData.coeff = squeeze(xSigMat(tt,rr,:,:))';
            xSigBoot{tt,rr} = getBootstrp(bootSigData,0,0.05);
            xSigBootAve(tt,rr,:) = xSigBoot{tt,rr}.coeff_bootave;
            % calculate two-tail p value
            for mm = 1:size(xCorrMat,3)

                if xCorrBoot{tt,rr}.coeff_bootave(mm) > 0
                    xCorrBootPval(tt,rr,mm) = sum(xCorrBoot{tt,rr}.bootSig(1,mm,:)<0)*2/size(xCorrBoot{tt,rr}.bootSig,3);
                else
                    xCorrBootPval(tt,rr,mm) = sum(xCorrBoot{tt,rr}.bootSig(1,mm,:)>0)*2/size(xCorrBoot{tt,rr}.bootSig,3);
                end

                 if xSigBoot{tt,rr}.coeff_bootave(mm) > 0
                    xSigBootPval(tt,rr,mm) = sum(xSigBoot{tt,rr}.bootSig(1,mm,:)<0)*2/size(xSigBoot{tt,rr}.bootSig,3);
                else
                    xSigBootPval(tt,rr,mm) = sum(xSigBoot{tt,rr}.bootSig(1,mm,:)>0)*2/size(xSigBoot{tt,rr}.bootSig,3);
                 end

            end

        end
    end
end

% plot crosscorrelation for every time bin
% add a star if significant from zero
colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    colorRange = [-0.1 0.1];
    t= -2.95:0.1:4.95;
 savefluofigpath = fullfile(save_path_fluo,'xcorrTime');
 if ~exist(savefluofigpath)
     mkdir(savefluofigpath);
 end

 varName = {'cn+1','cn','cn-1','rn','rn-1','xn','xn-1','ave_r','Cum_r','dQ','chosenQ','dK','chosenK','RPE','CKE'};


for pp = 1:size(xCorrMat,3)
     h=figure;

     subplot(1,2,1)
     image(xCorrBootAve(:,:,pp),'CDataMapping','scaled');
     axis square;
    colormap(colors);
    caxis([colorRange(1) colorRange(2)]);
    xticks(1:15)
    xticklabels(varName);
    xtickangle(45)
    yticks(1:15)
    yticklabels(varName)
    titleText = ['Cross-correlation at t= ',num2str(t(pp)),' s'];
    title(titleText);
    % add star for significance
    ifSig = xCorrBootPval(:,:,pp)<0.05;
     [r c]=find(ifSig==1);
    hold on;scatter(c,r,60,'MarkerEdgeColor',[254 177 57]/255,'marker','*');
    ax = gca(h);
    ax.FontSize = 20;
    %s.Marker = '*';
    subplot(3,30,90);
    image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
    colormap(colors);
    caxis([colorRange(1) colorRange(2)]);

    print(gcf,'-dpng',fullfile(savefluofigpath,['Cross-correlation_t=',num2str(pp)]));
    saveas(gcf, fullfile(savefluofigpath,['Cross-correlation_t=',num2str(pp)]), 'fig');
    saveas(gcf, fullfile(savefluofigpath,['Cross-correlation_t=',num2str(pp)]), 'svg');
    close;
end



% plot time curve as well

% get the average correlation first (from 0-3s)
xCorrBoot_time = cell(size(xCorrMat,1),size(xCorrMat,1));
xCorrBootAve_time =  zeros(size(xCorrMat,1),size(xCorrMat,1));
xCorrBootPval_time = ones(size(xCorrMat,1),size(xCorrMat,1));

xSigBoot_time = cell(size(xCorrMat,1),size(xCorrMat,1));
xSigBootAve_time =  zeros(size(xCorrMat,1),size(xCorrMat,1));
xSigBootPval_time = ones(size(xCorrMat,1),size(xCorrMat,1));

for tt = 1:size(xCorrMat,1)
    for rr = tt:size(xCorrMat,2)
        if tt~= rr
            bootData.coeff = squeeze(nanmean(xCorrMat(tt,rr,:,:),3));
            xCorrBoot_time{tt,rr} = getBootstrp(bootData,0,0.05);
            xCorrBootAve_time(tt,rr) = xCorrBoot_time{tt,rr}.coeff_bootave;

              bootSigData.coeff = squeeze(nanmean(xSigMat(tt,rr,:,:),3));
            xSigBoot_time{tt,rr} = getBootstrp(bootSigData,0,0.05);
            xSigBootAve_time(tt,rr) = xSigBoot_time{tt,rr}.coeff_bootave;
            % calculate two-tail p value


              bootSigData.coeff = squeeze(nanmean(xSigMat(tt,rr,:,:),3));
            xSigBoot_time{tt,rr} = getBootstrp(bootSigData,0,0.05);
            xSigBootAve_time(tt,rr) = xSigBoot_time{tt,rr}.coeff_bootave;

            if xCorrBoot_time{tt,rr}.coeff_bootave > 0
                xCorrBootPval_time(tt,rr) = sum(xCorrBoot_time{tt,rr}.bootSig(1,1,:)<0)*2/size(xCorrBoot_time{tt,rr}.bootSig,3);
            else
                xCorrBootPval_time(tt,rr) = sum(xCorrBoot_time{tt,rr}.bootSig(1,1,:)>0)*2/size(xCorrBoot_time{tt,rr}.bootSig,3);
            end

             if xSigBoot_time{tt,rr}.coeff_bootave > 0
                xSigBootPval_time(tt,rr) = sum(xSigBoot_time{tt,rr}.bootSig(1,1,:)<0)*2/size(xSigBoot_time{tt,rr}.bootSig,3);
            else
                xSigBootPval_time(tt,rr) = sum(xSigBoot_time{tt,rr}.bootSig(1,1,:)>0)*2/size(xSigBoot_time{tt,rr}.bootSig,3);
             end

        end

    end
end



% use bootstrap to check the significancy
varName = {'cn+1','cn','cn-1','rn','rn-1','xn','xn-1','ave_r','Cum_r','dQ','chosenQ','dK','chosenK','RPE','CKE'};


colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    colorRange = [-0.1 0.1];

    % coefficient xcorr
   h= figure;
    subplot(1,2,1)
   image(xCorrBootAve_time,'CDataMapping','scaled');
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
xticks(1:15)
xticklabels(varName)
xtickangle(45)
yticks(1:15)
yticklabels(varName)
ifSig = xCorrBootPval_time<0.05;
     [r c]=find(ifSig==1);
    hold on;scatter(c,r,60,'MarkerEdgeColor',[254 177 57]/255,'marker','*');
  ax = gca(h);
    ax.FontSize = 20;
    title('Average coeff cross-correlation');
    %s.Marker = '*';
    subplot(3,30,90);
    image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
    colormap(colors);
    caxis([colorRange(1) colorRange(2)]);


   print(gcf,'-dpng',fullfile(save_path_fluo,'Average coeff cross-correlation'));
    saveas(gcf, fullfile(save_path_fluo,'Average coeff cross-correlation'), 'fig');
    saveas(gcf, fullfile(save_path_fluo,'Average coeff cross-correlation'), 'svg');

     % significancy xcorr
   h= figure;
    subplot(1,2,1)
   image(xSigBootAve_time,'CDataMapping','scaled');
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
xticks(1:15)
xticklabels(varName)
xtickangle(45)
yticks(1:15)
yticklabels(varName)
ifSig = xSigBootPval_time<0.05;
     [r c]=find(ifSig==1);
    hold on;scatter(c,r,60,'MarkerEdgeColor',[254 177 57]/255,'marker','*');
  ax = gca(h);
    ax.FontSize = 20;
    title('Average sig cross-correlation');
    %s.Marker = '*';
    subplot(3,30,90);
    image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
    colormap(colors);
    caxis([colorRange(1) colorRange(2)]);



   print(gcf,'-dpng',fullfile(save_path_fluo,'Average cross-correlation'));
    saveas(gcf, fullfile(save_path_fluo,'Average cross-correlation'), 'fig');
    saveas(gcf, fullfile(save_path_fluo,'Average cross-correlation'), 'svg');

         % significancy xcorr
   h= figure;
    subplot(1,2,1)
   image(xSigBootAve_time,'CDataMapping','scaled');
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
xticks(1:15)
xticklabels(varName)
xtickangle(45)
yticks(1:15)
yticklabels(varName)
ifSig = xSigBootPval_time<0.05;
     [r c]=find(ifSig==1);
    hold on;scatter(c,r,60,'MarkerEdgeColor',[254 177 57]/255,'marker','*');
  ax = gca(h);
    ax.FontSize = 20;
    title('Average sig cross-correlation');
    %s.Marker = '*';
    subplot(3,30,90);
    image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
    colormap(colors);
    caxis([colorRange(1) colorRange(2)]);


   print(gcf,'-dpng',fullfile(save_path_fluo,'Average sig cross-correlation'));
    saveas(gcf, fullfile(save_path_fluo,'Average sig cross-correlation'), 'fig');
    saveas(gcf, fullfile(save_path_fluo,'Average sig cross-correlation'), 'svg');

%% negative crosscorrelation:
% possibility 1: different region, same modulation direction;
% possibility 2: same region, different direction; % check overlapping area

% negative pairs of interst:
% (dK/chosenK - dQ/chosenQ); cum_r;
% rn/RPE - ave_r;
% xn - cn/rn


% positive pairs of interst?
% dK(12) - cn+1(1)/cn-1(3); ave_r(8) (but not cn);

% other pairs:
% RPE(14)-cumR(9): positive -> negative
% rn(4) - dQ(10)/chosenK(13): negative for a period of time ~0.2-2s

% plot these curves
%% get the regression coefficient
RegCoeff = cell(1,15);
savematpath = fullfile(save_path_fluo,'Result');
reg1 = load(fullfile(savematpath,'regression1Sum.mat'));
%cn+1(2,cn(3,cn-1(4,rn(7,rn-1(8,xn(11.xn-1(12,ave_r(14,cum_r(15
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,2);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,2);coeff.boothigh=reg1.reg_cr_all.boothigh(:,2);
RegCoeff{1}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,3);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,3);coeff.boothigh=reg1.reg_cr_all.boothigh(:,3);
RegCoeff{2}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,4);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,4);coeff.boothigh=reg1.reg_cr_all.boothigh(:,4);
RegCoeff{3}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,7);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,7);coeff.boothigh=reg1.reg_cr_all.boothigh(:,7);
RegCoeff{4}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,8);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,8);coeff.boothigh=reg1.reg_cr_all.boothigh(:,8);
RegCoeff{5}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,11);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,11);coeff.boothigh=reg1.reg_cr_all.boothigh(:,11);
RegCoeff{6}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,12);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,12);coeff.boothigh=reg1.reg_cr_all.boothigh(:,12);
RegCoeff{7}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,14);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,14);coeff.boothigh=reg1.reg_cr_all.boothigh(:,14);
RegCoeff{8}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,15);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,15);coeff.boothigh=reg1.reg_cr_all.boothigh(:,15);
RegCoeff{9}=coeff;
%dQ(8,chosenQ(9,dK(10,chosenK(11
reg1 = load(fullfile(savematpath,'regression2Sum.mat'));
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,8);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,8);coeff.boothigh=reg1.reg_cr_all.boothigh(:,8);
RegCoeff{10}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,9);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,9);coeff.boothigh=reg1.reg_cr_all.boothigh(:,9);
RegCoeff{11}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,10);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,10);coeff.boothigh=reg1.reg_cr_all.boothigh(:,10);
RegCoeff{12}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,11);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,11);coeff.boothigh=reg1.reg_cr_all.boothigh(:,11);
RegCoeff{13}=coeff;
%RPE(6,CKE(8
reg1 = load(fullfile(savematpath,'regression3Sum.mat'));
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,6);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,6);coeff.boothigh=reg1.reg_cr_all.boothigh(:,6);
RegCoeff{14}=coeff;
coeff.coeff_bootave=reg1.reg_cr_all.coeff_bootave(:,8);
coeff.bootlow=reg1.reg_cr_all.bootlow(:,8);coeff.boothigh=reg1.reg_cr_all.boothigh(:,8);
RegCoeff{15}=coeff;

%% plot
% RPE & CumR
t=-2.95:0.1:4.95;
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,9,14,varName,save_path_fluo);

% rn & dQ
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,4,10,varName,save_path_fluo);
% rn & chosenK
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,4,13,varName,save_path_fluo);
% dK & cn-1
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,3,12,varName,save_path_fluo);
% dK & cn+1
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,1,12,varName,save_path_fluo);
% dK & cn
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,2,12,varName,save_path_fluo);
% dK & ave_r
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,8,12,varName,save_path_fluo);
% rn & ave_r
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,4,8,varName,save_path_fluo);
% rn & cum_r
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,4,9,varName,save_path_fluo);
% cn & rn
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,2,4,varName,save_path_fluo);


% rn & xn
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,4,6,varName,save_path_fluo);
% cn & xn
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,2,6,varName,save_path_fluo);
% dK & dQ
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,10,12,varName,save_path_fluo);
% dK & cumR
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,9,12,varName,save_path_fluo);

% rn & RPE
plot_pairedXcorr(xCorrBoot,RegCoeff,sigVar,t,4,14,varName,save_path_fluo);

% dQ & cumR
%% check the result: calculate overlap area? - plot average coefficient
% negative crosscorrelation:
% possibility 1: different region, same modulation direction;
% possibility 2: same region, different direction; % check overlapping area

end
