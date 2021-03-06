function MP_GRAB_temporalCorrSummary(dataIndex, savesumfigpath,savematsumpath)

% calculate temporal correlation of different varible coefficients within
% grids

% check if coefficient sign changes over time
% for eachgrid, calculate the number of + time point and - time point
% then average

nFiles = size(dataIndex,1);

%% get clustering results, average similar clusters, also rerun clustering on all grids together
rCoeff = cell(0); rClust = cell(0);

% combine isSig list, for venn plot
SigList.c= [];SigList.r=[];SigList.x=[];
SigList.dQ=[];SigList.dK=[];SigList.pRPE=[];SigList.nRPE=[];
SigList.CKE=[];
SigList.session=[];
% count the number of significant grids for every variables and sessions
numSigVar = zeros(8,nFiles);
numSigTotal = zeros(1,nFiles);
varList = {'choice', 'outcome', 'interaction', 'dQ', 'dK', 'posRPE', 'negRPE', 'CKE'};
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
        
        % LOAD SAVED DATA
        saveregpath = fullfile(savematpath,'cluster.mat');
        if exist(saveregpath)
            load(saveregpath)
        end
%         % outcome
%         tlabel1 = 'Outcome coefficient';
%         saveregpath =  fullfile(savematpath,[tlabel1,' cluster.mat']);  % file path to save the results
%         if exist(saveregpath)
%             load(saveregpath)
%         end

        rCoeff{ii} = saveDataOutcome.coeff(saveDataOutcome.oriInd,:);
        rClust{ii} = saveDataOutcome.T;
        rOriID{ii} = saveDataOutcome.oriInd;
        numSigVar(2,ii) = length(saveDataOutcome.oriInd);
        rt = saveDataOutcome.t;
        
        choiceCoeff{ii} = choicetempData.coeff;choiceisSig{ii} = choicetempData.isSig; numSigVar(1,ii) = sum(choicetempData.isSig);
        xnCoeff{ii} = xntempData.coeff;xnisSig{ii} = xntempData.isSig;numSigVar(3,ii) = sum(xntempData.isSig);
        posRPECoeff{ii} = posRPEtempData.coeff; posRPEisSig{ii} = posRPEtempData.isSig;numSigVar(6,ii) = sum(posRPEtempData.isSig);
        negRPECoeff{ii} = negRPEtempData.coeff; negRPEisSig{ii} = negRPEtempData.isSig;numSigVar(7,ii) = sum(posRPEtempData.isSig);
        %negRPECoeff{ii} = negRPEtempData.coeff; negRPEisSig{ii} = negRPEtempData.isSig;
        dQCoeff{ii} = dQtempData.coeff; dQisSig{ii} = dQtempData.isSig;numSigVar(4,ii) = sum(dQtempData.isSig);
        dKCoeff{ii} = dKtempData.coeff; dKisSig{ii} = dKtempData.isSig;numSigVar(5,ii) = sum(dKtempData.isSig);
        CKECoeff{ii} = CKEtempData.coeff; CKEisSig{ii} = CKEtempData.isSig;numSigVar(8,ii) = sum(CKEtempData.isSig);
        
        % sig list for venn plot
        SigList.session = [SigList.session,ii*ones(1,length(choicetempData.isSig))];
        SigList.c= [SigList.c,choicetempData.isSig'];
        outcomeSig = zeros(1,size(saveDataOutcome.coeff,1));
        outcomeSig(saveDataOutcome.sigInd)=1;
        SigList.r = [SigList.r, outcomeSig];
        SigList.x = [SigList.x, xntempData.isSig'];
        SigList.dQ = [SigList.dQ, dQtempData.isSig'];
        SigList.dK = [SigList.dK, dKtempData.isSig'];
        SigList.pRPE = [SigList.pRPE, posRPEtempData.isSig'];
        SigList.nRPE = [SigList.nRPE, negRPEtempData.isSig'];
        SigList.CKE = [SigList.CKE, CKEtempData.isSig'];

        numSigTotal(ii) = sum(choicetempData.isSig|outcomeSig'|xntempData.isSig|dQtempData.isSig|dKtempData.isSig|posRPEtempData.isSig|negRPEtempData.isSig|CKEtempData.isSig);
    end
end

% save the siglist
if ~exist(savematsumpath)
    mkdir(savematsumpath);
end
save(fullfile(savematsumpath,'sigList.mat'),'SigList');
save(fullfile(savematsumpath,'sigNum.mat'),'numSigVar','numSigTotal');

%% plot percentage of significant grids
sigForAny = SigList.c|SigList.r|SigList.x|SigList.dQ|SigList.dK|SigList.pRPE|SigList.nRPE|SigList.CKE;
numSigAny = sum(sigForAny);
figure; 
explode = [0 1];
p=pie([numSigAny,length(sigForAny)-numSigAny],explode);
print(gcf,'-dpng',fullfile(savesumfigpath,'Percentage of significant grids All'));
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids All'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids All'), 'svg');

% violin plot for significant grids with different variables
figure;violinplot(numSigVar'./size(choiceCoeff{1},1),varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
ylabel('Fraction of significant grids');
print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables'));
saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables'), 'svg');

% 

%% venn diagram showing overlap 

allInd = 1:length(SigList.c);
cSigInd = allInd(logical(SigList.c)); rSigInd = allInd(logical(SigList.r)); xSigInd = allInd(logical(SigList.x));
dQSigInd = allInd(logical(SigList.dQ)); dKSigInd = allInd(logical(SigList.dK));
pRPESigInd = allInd(logical(SigList.pRPE)); nRPESigInd = allInd(logical(SigList.nRPE));
CKESigInd = allInd(logical(SigList.CKE));

%% test for overlapping significance

% figure;
% h = vennEulerDiagram({dQSigInd;dKSigInd;pRPESigInd;nRPESigInd}, 'drawProportional', true, 'SetLabels', ["dQ"; "dK"; "posRPE"; "negRPE"]);
% 
% figure;
% h = vennEulerDiagram({rSigInd;dQSigInd;pRPESigInd;nRPESigInd}, 'drawProportional', true, 'SetLabels', ["Outcome"; "dQ";"posRPE"; "negRPE"]);
% 
% figure;
% h = vennEulerDiagram({dQSigInd;pRPESigInd;nRPESigInd}, 'drawProportional', true, 'SetLabels', ["dQ"; "posRPE"; "negRPE"]);
% 
% figure;
% h = vennEulerDiagram({cSigInd;xSigInd;dKSigInd}, 'drawProportional', true, 'SetLabels', ["choice"; "interaction"; "dK"]);
% 
% figure;
% h = vennEulerDiagram({cSigInd;xSigInd;dQSigInd;dKSigInd}, 'drawProportional', true, 'SetLabels', ["choice"; "interaction"; 'dQ'; 'dK']);

% chi-square test for independence

% choice and outcome
x = [length(intersect(cSigInd,rSigInd)),length(setdiff(cSigInd,rSigInd));
    length(setdiff(rSigInd,cSigInd)),length(allInd)-length(union(cSigInd,rSigInd))];
[h,chi,p]=chi2ind(x,0.05)
y=[length(intersect(cSigInd,xSigInd)),length(setdiff(cSigInd,xSigInd));
    length(setdiff(xSigInd,cSigInd)),length(allInd)-length(union(cSigInd,xSigInd))];
[h,chi,p]=chi2ind(y,0.05)
z=[length(intersect(rSigInd,xSigInd)),length(setdiff(rSigInd,xSigInd));
    length(setdiff(xSigInd,rSigInd)),length(allInd)-length(union(rSigInd,xSigInd))];
[h,chi,p]=chi2ind(z,0.05)
close all;
%% trim the clustering.
% 1. delete cluster with too few grids . <10 grids?
% 2. combine goup with high correlation >0.9?
g1=cell(0); g2 = cell(0); g3 = cell(0);
clustCoeff = []; clustMean = [];
rClust_new = cell(0);
rCoeff_new = cell(0);
cCoeff_new = cell(0);
xCoeff_new = cell(0);
posRPECoeff_new = cell(0);
negRPECoeff_new = cell(0);
dQCoeff_new = cell(0);
dKCoeff_new = cell(0);
CKECoeff_new = cell(0);

% only take significant cells
% cCoeff_Sig = cell(0);
% xCoeff_Sig = cell(0);
% posRPECoeff_Sig = cell(0);
% negRPECoeff_Sig = cell(0);

% sig ROIs that not sig for outcome
cCoeff_nR = [];
xCoeff_nR = [];
posRPECoeff_nR = [];
negRPECoeff_nR = [];
dQCoeff_nR = [];
dKCoeff_nR = [];
CKECoeff_nR = [];

rOriID_new = cell(0);
newInd = 1;
% rHPeakT = [];
% rLPeakT = [];
% rDecayTau = [];
oriGroup = [];
% save all grids that cannot be grouped into the 3 groups



t0 = 0; t1 = 4;  % time range for correlation calculation
ROIId = 1:size(choiceCoeff{1},1);
% for sig grids counting for group4 (not sig for outcome)
numSigVarg4 = zeros(8,nFiles);

for gg = 1:nFiles
    g1Num = sum(rClust{gg}==1); g2Num = sum(rClust{gg}==2); g3Num = sum(rClust{gg}==3);
    % calculate the average coefficient with in 1-2.5s for each group
    g1Mean =nanmean(rCoeff{gg}(rClust{gg}==1,:),1);  g2Mean = nanmean(rCoeff{gg}(rClust{gg}==2,:),1); g3Mean =nanmean(rCoeff{gg}(rClust{gg}==3,:),1);

    %clustMean = [clustMean;g1Mean;g2Mean;g3Mean];
    %% smooth the mean and find peak and decay tau
    % find 1): highest peak and time;
    % 2): lowest peak and time;
    % 3) exponential decay tau after high peak



    % (islocalmin & islocalmax)
    tInd = 1:length(rt);
    searchInd = tInd(rt>0);
    % group 1
%     smoothed1 = smooth(g1Mean);
%     [minValue1,minInd] = nanmin(smoothed1(rt>0));
%     minInd1 = searchInd(minInd);
%     [maxValue1,maxInd] = nanmax(smoothed1(rt>0));
%     maxInd1 = searchInd(maxInd);
%     rHPeakT = [rHPeakT, rt(maxInd1)];
%     rLPeakT = [rLPeakT, rt(minInd1)];
    oriGroup = [oriGroup,1];
    % fit for tau
    %     if maxInd1<length(rt)
    %         f1 = fit(rt(maxInd1:end),smoothed1(maxInd1:end),'exp1');
    %         rDecayTau = [rDecayTau, f1.b];
    %     else
    %         f1 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==1,:);
    cCoeff_new{newInd} = choiceCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    xCoeff_new{newInd} = xnCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    posRPECoeff_new{newInd} = posRPECoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    negRPECoeff_new{newInd} = negRPECoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    dQCoeff_new{newInd} = dQCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    dKCoeff_new{newInd} = dKCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    CKECoeff_new{newInd} = CKECoeff{gg}(rOriID{gg}(rClust{gg}==1),:);

    group1Sig = zeros(size(choiceCoeff{1},1),1);
    group1SigId = rOriID{gg}(rClust{gg}==1);
    choiceSigId = ROIId(choiceisSig{gg});
    xnSigId = ROIId(xnisSig{gg});
    posRPESigId = ROIId(posRPEisSig{gg});
    negRPESigId = ROIId(negRPEisSig{gg});
    dQSigId = ROIId(dQisSig{gg});
    dKSigId = ROIId(dKisSig{gg});
    CKESigId = ROIId(CKEisSig{gg});

    choiceSig = zeros(1,length(group1SigId));
    xnSig = zeros(1,length(group1SigId));
    posRPESig = zeros(1,length(group1SigId));
    negRPESig = zeros(1,length(group1SigId));
    dQSig = zeros(1,length(group1SigId));
    dKSig = zeros(1,length(group1SigId));
    CKESig = zeros(1,length(group1SigId));

    for rr = 1:length(group1SigId)
        if ismember(group1SigId(rr),choiceSigId)
            choiceSig(rr) = 1;
        end
        if ismember(group1SigId(rr),xnSigId)
            xnSig(rr) = 1;
        end
        if ismember(group1SigId(rr),posRPESigId)
            posRPESig(rr) = 1;
        end
        if ismember(group1SigId(rr),negRPESigId)
            negRPESig(rr) = 1;
        end
        if ismember(group1SigId(rr),dQSigId)
            dQSig(rr) = 1;
        end
        if ismember(group1SigId(rr),dKSigId)
            dKSig(rr) = 1;
        end
        if ismember(group1SigId(rr),CKESigId)
            CKESig(rr) = 1;
        end
    end
    cisSig{newInd} = choiceSig;
    xisSig{newInd} = xnSig;
    pRPEisSig{newInd} = posRPESig;
    nRPEisSig{newInd} = negRPESig;
    dqisSig{newInd} = dQSig;
    dkisSig{newInd} = dKSig;
    ckeisSig{newInd} = CKESig;

%     cCoeff_Sig{newInd} = choiceCoeff{gg}(group1Sig & choiceisSig{gg},:);
%     xCoeff_Sig{newInd} = xnCoeff{gg}(group1Sig & xnisSig{gg},:);
%     posRPECoeff_Sig{newInd} = posRPECoeff{gg}(group1Sig & posRPEisSig{gg},:);
%     negRPECoeff_Sig{newInd} = negRPECoeff{gg}(group1Sig & negRPEisSig{gg},:);
%     
    RisSig = zeros(1,size(choiceCoeff{1},1));
    RisSig(rOriID{gg}) = 1;

    cCoeff_nR = [cCoeff_nR; choiceCoeff{gg}(~RisSig'&choiceisSig{gg},:)];numSigVarg4(1,gg) = sum(~RisSig'&choiceisSig{gg});
    xCoeff_nR = [xCoeff_nR; xnCoeff{gg}(~RisSig'&xnisSig{gg},:)];numSigVarg4(3,gg) = sum(~RisSig'&xnisSig{gg});
    posRPECoeff_nR = [posRPECoeff_nR; posRPECoeff{gg}(~RisSig'&posRPEisSig{gg},:)];numSigVarg4(6,gg) = sum(~RisSig'&posRPEisSig{gg});
    negRPECoeff_nR = [negRPECoeff_nR; negRPECoeff{gg}(~RisSig'&negRPEisSig{gg},:)];numSigVarg4(7,gg) = sum(~RisSig'&negRPEisSig{gg});
    dQCoeff_nR = [dQCoeff_nR; dQCoeff{gg}(~RisSig'&dQisSig{gg},:)];numSigVarg4(4,gg) = sum(~RisSig'&dQisSig{gg});
    dKCoeff_nR = [dKCoeff_nR; dKCoeff{gg}(~RisSig'&dKisSig{gg},:)];numSigVarg4(5,gg) = sum(~RisSig'&dKisSig{gg});
    CKECoeff_nR = [CKECoeff_nR; CKECoeff{gg}(~RisSig'&CKEisSig{gg},:)];numSigVarg4(8,gg) = sum(~RisSig'&CKEisSig{gg});

    newInd = newInd+1;

    % group 2
%     smoothed2 = smooth(g2Mean);
%     [minValue2,minInd] = nanmin(smoothed2(rt>0));
%     minInd2 = searchInd(minInd);
%     [maxValue2,maxInd] = nanmax(smoothed2(rt>0));
%     maxInd2 = searchInd(maxInd);
%     rHPeakT = [rHPeakT, rt(maxInd2)];
%     rLPeakT = [rLPeakT, rt(minInd2)];
    oriGroup = [oriGroup,2];
    % fit for tau
    %     if maxInd2<length(rt)
    %         f2 = fit(rt(maxInd2:end),smoothed2(maxInd2:end),'exp1');
    %         rDecayTau = [rDecayTau, f2.b];
    %     else
    %         f2 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==2,:);
     cCoeff_new{newInd} = choiceCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    xCoeff_new{newInd} = xnCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    posRPECoeff_new{newInd} = posRPECoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    negRPECoeff_new{newInd} = negRPECoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    dQCoeff_new{newInd} = dQCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    dKCoeff_new{newInd} = dKCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    CKECoeff_new{newInd} = CKECoeff{gg}(rOriID{gg}(rClust{gg}==2),:);

    group2Sig = zeros(size(choiceCoeff{1},1),1);
    group2Sig(rOriID{gg}(rClust{gg}==2)) = 1;
    group2SigId = rOriID{gg}(rClust{gg}==2);

    choiceSig = zeros(1,length(group2SigId));
    xnSig = zeros(1,length(group2SigId));
    posRPESig = zeros(1,length(group2SigId));
    negRPESig = zeros(1,length(group2SigId));
    dQSig = zeros(1,length(group2SigId));
    dKSig = zeros(1,length(group2SigId));
    CKESig = zeros(1,length(group2SigId));

    for rr = 1:length(group2SigId)
        if ismember(group2SigId(rr),choiceSigId)
            choiceSig(rr) = 1;
        end
        if ismember(group2SigId(rr),xnSigId)
            xnSig(rr) = 1;
        end
        if ismember(group2SigId(rr),posRPESigId)
            posRPESig(rr) = 1;
        end
        if ismember(group2SigId(rr),negRPESigId)
            negRPESig(rr) = 1;
        end
        if ismember(group2SigId(rr),dQSigId)
            dQSig(rr) = 1;
        end
        if ismember(group2SigId(rr),dKSigId)
            dKSig(rr) = 1;
        end
        if ismember(group2SigId(rr),CKESigId)
            CKESig(rr) = 1;
        end
    end
    cisSig{newInd} = choiceSig;
    xisSig{newInd} = xnSig;
    pRPEisSig{newInd} = posRPESig;
    nRPEisSig{newInd} = negRPESig;
    dqisSig{newInd} = dQSig;
    dkisSig{newInd} = dKSig;
    ckeisSig{newInd} = CKESig;
%      cCoeff_Sig{newInd} = choiceCoeff{gg}(group2Sig & choiceisSig{gg},:);
%      xCoeff_Sig{newInd} = xnCoeff{gg}(group2Sig & xnisSig{gg},:);
%     posRPECoeff_Sig{newInd} = posRPECoeff{gg}(group2Sig & posRPEisSig{gg},:);
%     negRPECoeff_Sig{newInd} = negRPECoeff{gg}(group2Sig & negRPEisSig{gg},:);
%     
    newInd = newInd+1;

    % group 3
%     smoothed3 = smooth(g3Mean);
%     [minValue3,minInd] = nanmin(smoothed3(rt>0));
%     minInd3 = searchInd(minInd);
%     [maxValue3,maxInd] = nanmax(smoothed3(rt>0));
%     maxInd3 = searchInd(maxInd);
%     rHPeakT = [rHPeakT, rt(maxInd3)];
%     rLPeakT = [rLPeakT, rt(minInd3)];
    oriGroup = [oriGroup,3];
    % fit for tau
    %     if maxInd3<length(rt)
    %         f3 = fit(rt(maxInd3:end),smoothed3(maxInd3:end),'exp1');
    %         rDecayTau = [rDecayTau, f3.b];
    %     else
    %         f3 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==3,:);
     cCoeff_new{newInd} = choiceCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    xCoeff_new{newInd} = xnCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    posRPECoeff_new{newInd} = posRPECoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    negRPECoeff_new{newInd} = negRPECoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    dQCoeff_new{newInd} = dQCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    dKCoeff_new{newInd} = dKCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    CKECoeff_new{newInd} = CKECoeff{gg}(rOriID{gg}(rClust{gg}==3),:);

     group3Sig = zeros(size(choiceCoeff{1},1),1);
    group3Sig(rOriID{gg}(rClust{gg}==3)) = 1;
     group3SigId = rOriID{gg}(rClust{gg}==3);

    choiceSig = zeros(1,length(group3SigId));
    xnSig = zeros(1,length(group3SigId));
    posRPESig = zeros(1,length(group3SigId));
    negRPESig = zeros(1,length(group3SigId));
    dQSig = zeros(1,length(group3SigId));
    dKSig = zeros(1,length(group3SigId));
    CKESig = zeros(1,length(group3SigId));

    for rr = 1:length(group3SigId)
        if ismember(group3SigId(rr),choiceSigId)
            choiceSig(rr) = 1;
        end
        if ismember(group3SigId(rr),xnSigId)
            xnSig(rr) = 1;
        end
        if ismember(group3SigId(rr),posRPESigId)
            posRPESig(rr) = 1;
        end
        if ismember(group3SigId(rr),negRPESigId)
            negRPESig(rr) = 1;
        end
        if ismember(group3SigId(rr),dQSigId)
            dQSig(rr) = 1;
        end
        if ismember(group3SigId(rr),dKSigId)
            dKSig(rr) = 1;
        end
        if ismember(group3SigId(rr),CKESigId)
            CKESig(rr) = 1;
        end
    end
    cisSig{newInd} = choiceSig;
    xisSig{newInd} = xnSig;
    pRPEisSig{newInd} = posRPESig;
    nRPEisSig{newInd} = negRPESig;
    dqisSig{newInd} = dQSig;
    dkisSig{newInd} = dKSig;
    ckeisSig{newInd} = CKESig;

%      cCoeff_Sig{newInd} = choiceCoeff{gg}(group3Sig & choiceisSig{gg},:);
%      xCoeff_Sig{newInd} = xnCoeff{gg}(group3Sig & xnisSig{gg},:);
%      posRPECoeff_Sig{newInd} = posRPECoeff{gg}(group3Sig & posRPEisSig{gg},:);
%      negRPECoeff_Sig{newInd} = negRPECoeff{gg}(group3Sig & negRPEisSig{gg},:);
%     
    newInd = newInd+1;

    % plot the fit results
    savefigpath =  fullfile(dataIndex.BehPath{gg},'figs-fluo','clustProp');
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end
  

end

% regroup with criteria *area under curve
group1 = []; group2 = []; group3 = []; %group4 = [];
choice1 = []; choice2 = []; choice3 = [];
cSig1 = []; cSig2 = []; cSig3 = [];
%choice1Sig = [];choice2Sig = [];choice3Sig = [];
xn1 = []; xn2 = []; xn3 = [];
xSig1 = []; xSig2 = []; xSig3 = [];
%xn1Sig = []; xn2Sig = []; xn3Sig = []; 
posRPE1 = []; posRPE2 = []; posRPE3 = [];
pRPESig1 = []; pRPESig2 = []; pRPESig3 = [];
%posRPE1Sig = []; posRPE2Sig = []; posRPE3Sig = []; 
negRPE1 = []; negRPE2 = []; negRPE3 = [];
nRPESig1 = []; nRPESig2 = []; nRPESig3 = [];
%negRPE1Sig = []; negRPE2Sig = []; negRPE3Sig = []; 
dQ1 = []; dQ2 = []; dQ3 = [];
dQSig1 = []; dQSig2 = []; dQSig3 = [];

dK1 = []; dK2 = []; dK3 = [];
dKSig1 = []; dKSig2 = []; dKSig3 = [];

CKE1 = []; CKE2 = []; CKE3 = [];
CKESig1 = []; CKESig2 = []; CKESig3 = [];
% if the surface area between (0,2) less than 0, group1
% if the surface area between (2,4) less than 0, group4

%
numSigVarg1 = zeros(8,nFiles);
numSigVarg2 = zeros(8,nFiles);
numSigVarg3 = zeros(8,nFiles);

for ggg = 1:length(rCoeff_new)
    animalInd = floor((ggg-1)/3)+1;
    if oriGroup(ggg) == 1
        line11 = nanmean(rCoeff_new{ggg},1);
        areaUnder1 = trapz(rt(rt>0&rt<2),line11(rt>0 & rt<2));
       %areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
        if areaUnder1 < 0 %& areaUnder2 > areaUnder1
            group1 = [group1;rCoeff_new{ggg}];
            numSigVarg1(2,animalInd) = numSigVarg1(2,animalInd) + size(rCoeff_new{ggg},1);
            choice1 = [choice1; cCoeff_new{ggg}]; cSig1 = [cSig1,cisSig{ggg}];%choice1Sig = [choice1Sig; cCoeff_Sig{ggg}];
            numSigVarg1(1,animalInd) = numSigVarg1(1,animalInd) + sum(cisSig{ggg});
            xn1 = [xn1; xCoeff_new{ggg}]; xSig1 = [xSig1,xisSig{ggg}];%xn1Sig = [xn1Sig; xCoeff_Sig{ggg}];
            numSigVarg1(3,animalInd) = numSigVarg1(3,animalInd) + sum(xisSig{ggg});
            posRPE1 = [posRPE1; posRPECoeff_new{ggg}]; pRPESig1 = [pRPESig1,pRPEisSig{ggg}]; %posRPE1Sig = [posRPE1Sig; posRPECoeff_Sig{ggg}];
            numSigVarg1(6,animalInd) = numSigVarg1(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE1 = [negRPE1; negRPECoeff_new{ggg}]; nRPESig1 = [nRPESig1,nRPEisSig{ggg}];%negRPE1Sig = [negRPE1Sig; negRPECoeff_Sig{ggg}];
            numSigVarg1(7,animalInd) = numSigVarg1(7,animalInd) + sum(nRPEisSig{ggg});
            dQ1 = [dQ1; dQCoeff_new{ggg}]; dQSig1 = [dQSig1,dqisSig{ggg}];
            numSigVarg1(4,animalInd) = numSigVarg1(4,animalInd) + sum(dqisSig{ggg});
            dK1 = [dK1; dKCoeff_new{ggg}]; dKSig1 = [dKSig1,dkisSig{ggg}];
            numSigVarg1(5,animalInd) = numSigVarg1(5,animalInd) + sum(dkisSig{ggg});
            CKE1 = [CKE1; CKECoeff_new{ggg}]; CKESig1 = [CKESig1,ckeisSig{ggg}];
            numSigVarg1(8,animalInd) = numSigVarg1(8,animalInd) + sum(ckeisSig{ggg});
            
        else
            group2 = [group2;rCoeff_new{ggg}];
            numSigVarg2(2,animalInd) = numSigVarg2(2,animalInd) + size(rCoeff_new{ggg},1);
            choice2 = [choice2; cCoeff_new{ggg}]; cSig2 = [cSig2,cisSig{ggg}];%choice2Sig = [choice2Sig; cCoeff_Sig{ggg}];
            numSigVarg2(1,animalInd) = numSigVarg2(1,animalInd) + sum(cisSig{ggg});
            xn2 = [xn2; xCoeff_new{ggg}]; xSig2 = [xSig2,xisSig{ggg}];%xn2Sig = [xn2Sig; xCoeff_Sig{ggg}];
            numSigVarg2(3,animalInd) = numSigVarg2(3,animalInd) + sum(xisSig{ggg});
            posRPE2 = [posRPE2; posRPECoeff_new{ggg}]; pRPESig2 = [pRPESig2,pRPEisSig{ggg}]; %posRPE2Sig = [posRPE2Sig; posRPECoeff_Sig{ggg}];
            numSigVarg2(6,animalInd) = numSigVarg2(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE2 = [negRPE2; negRPECoeff_new{ggg}]; nRPESig2 = [nRPESig2,nRPEisSig{ggg}];%negRPE2Sig = [negRPE2Sig; negRPECoeff_Sig{ggg}];
            numSigVarg2(7,animalInd) = numSigVarg2(7,animalInd) + sum(nRPEisSig{ggg}); 
            dQ2 = [dQ2; dQCoeff_new{ggg}]; dQSig2 = [dQSig2,dqisSig{ggg}];
            numSigVarg2(4,animalInd) = numSigVarg2(4,animalInd) + sum(dqisSig{ggg});
             dK2 = [dK2; dKCoeff_new{ggg}]; dKSig2 = [dKSig2,dkisSig{ggg}];
             numSigVarg2(5,animalInd) = numSigVarg2(5,animalInd) + sum(dkisSig{ggg});
            CKE2 = [CKE2; CKECoeff_new{ggg}]; CKESig2 = [CKESig2,ckeisSig{ggg}];
           numSigVarg2(8,animalInd) = numSigVarg2(8,animalInd) + sum(ckeisSig{ggg});
        end
    elseif oriGroup(ggg) == 2
        line11 = nanmean(rCoeff_new{ggg},1);
        areaUnder1 = trapz(rt(rt>0&rt<2),line11(rt>0 & rt<2));
        %areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
        if areaUnder1 < 0 %& areaUnder2 > areaUnder1
            group1 = [group1;rCoeff_new{ggg}];
            numSigVarg1(2,animalInd) = numSigVarg1(2,animalInd) + size(rCoeff_new{ggg},1);
%         elseif areaUnder2 < 0
%             group4 = [group4;rCoeff_new{ggg}];
            choice1 = [choice1; cCoeff_new{ggg}]; cSig1 = [cSig1,cisSig{ggg}];%choice1Sig = [choice1Sig; cCoeff_Sig{ggg}];
            numSigVarg1(1,animalInd) = numSigVarg1(1,animalInd) + sum(cisSig{ggg});
            xn1 = [xn1; xCoeff_new{ggg}]; xSig1 = [xSig1,xisSig{ggg}];%xn1Sig = [xn1Sig; xCoeff_Sig{ggg}];
            numSigVarg1(3,animalInd) = numSigVarg1(3,animalInd) + sum(xisSig{ggg});
            posRPE1 = [posRPE1; posRPECoeff_new{ggg}]; pRPESig1 = [pRPESig1,pRPEisSig{ggg}];%posRPE1Sig = [posRPE1Sig; posRPECoeff_Sig{ggg}];
            numSigVarg1(6,animalInd) = numSigVarg1(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE1 = [negRPE1; negRPECoeff_new{ggg}]; nRPESig1 = [nRPESig1,nRPEisSig{ggg}];%negRPE2Sig = [negRPE2Sig; negRPECoeff_Sig{ggg}];
            numSigVarg1(7,animalInd) = numSigVarg1(7,animalInd) + sum(nRPEisSig{ggg});  
            dQ1 = [dQ1; dQCoeff_new{ggg}]; dQSig1 = [dQSig1,dqisSig{ggg}];
             numSigVarg1(4,animalInd) = numSigVarg1(4,animalInd) + sum(dqisSig{ggg});
            dK1 = [dK1; dKCoeff_new{ggg}]; dKSig1 = [dKSig1,dkisSig{ggg}];
            numSigVarg1(5,animalInd) = numSigVarg1(5,animalInd) + sum(dkisSig{ggg});
            CKE1 = [CKE1; CKECoeff_new{ggg}]; CKESig1 = [CKESig1,ckeisSig{ggg}];
           numSigVarg1(8,animalInd) = numSigVarg1(8,animalInd) + sum(ckeisSig{ggg});
            
        else
            group2 = [group2;rCoeff_new{ggg}];
            numSigVarg2(2,animalInd) = numSigVarg2(2,animalInd) + size(rCoeff_new{ggg},1);
            choice2 = [choice2; cCoeff_new{ggg}]; cSig2 = [cSig2,cisSig{ggg}];%choice2Sig = [choice2Sig; cCoeff_Sig{ggg}];
            numSigVarg2(1,animalInd) = numSigVarg2(1,animalInd) + sum(cisSig{ggg});
            xn2 = [xn2; xCoeff_new{ggg}]; xSig2 = [xSig2,xisSig{ggg}];%xn2Sig = [xn2Sig; xCoeff_Sig{ggg}];
            numSigVarg2(3,animalInd) = numSigVarg2(3,animalInd) + sum(xisSig{ggg});
            posRPE2 = [posRPE2; posRPECoeff_new{ggg}]; pRPESig2 = [pRPESig2,pRPEisSig{ggg}];%posRPE2Sig = [posRPE2Sig; posRPECoeff_Sig{ggg}];
            numSigVarg2(6,animalInd) = numSigVarg2(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE2 = [negRPE2; negRPECoeff_new{ggg}]; nRPESig2 = [nRPESig2,nRPEisSig{ggg}];%negRPE2Sig = [negRPE2Sig; negRPECoeff_Sig{ggg}];
            numSigVarg2(7,animalInd) = numSigVarg2(7,animalInd) + sum(nRPEisSig{ggg});
            dQ2 = [dQ2; dQCoeff_new{ggg}]; dQSig2 = [dQSig2,dqisSig{ggg}];
            numSigVarg2(4,animalInd) = numSigVarg2(4,animalInd) + sum(dqisSig{ggg});
            dK2 = [dK2; dKCoeff_new{ggg}]; dKSig2 = [dKSig2,dkisSig{ggg}];
             numSigVarg2(5,animalInd) = numSigVarg2(5,animalInd) + sum(dkisSig{ggg});
            CKE2 = [CKE2; CKECoeff_new{ggg}]; CKESig2 = [CKESig2,ckeisSig{ggg}];
              numSigVarg2(8,animalInd) = numSigVarg2(8,animalInd) + sum(ckeisSig{ggg});
        end
    else
%         areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
%         if areaUnder2 < 0
%             group4 = [group4;rCoeff_new{ggg}];
%         else
            group3 = [group3;rCoeff_new{ggg}];
            numSigVarg3(2,animalInd) = numSigVarg3(2,animalInd) + size(rCoeff_new{ggg},1);
            choice3 = [choice3; cCoeff_new{ggg}]; cSig3 = [cSig3,cisSig{ggg}];
            numSigVarg3(1,animalInd) = numSigVarg3(1,animalInd) + sum(cisSig{ggg});%choice3Sig = [choice3Sig; cCoeff_Sig{ggg}];
            xn3 = [xn3; xCoeff_new{ggg}]; xSig3 = [xSig3,xisSig{ggg}];%xn3Sig = [xn3Sig; xCoeff_Sig{ggg}];
            numSigVarg3(3,animalInd) = numSigVarg3(3,animalInd) + sum(xisSig{ggg});
            posRPE3 = [posRPE3; posRPECoeff_new{ggg}]; pRPESig3 = [pRPESig3,pRPEisSig{ggg}];%posRPE3Sig = [posRPE3Sig; posRPECoeff_Sig{ggg}];
            numSigVarg3(6,animalInd) = numSigVarg3(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE3 = [negRPE3; negRPECoeff_new{ggg}]; nRPESig3 = [nRPESig3,nRPEisSig{ggg}];%negRPE3Sig = [negRPE3Sig; negRPECoeff_Sig{ggg}];
            numSigVarg3(7,animalInd) = numSigVarg3(7,animalInd) + sum(nRPEisSig{ggg});
            dQ3 = [dQ3; dQCoeff_new{ggg}]; dQSig3 = [dQSig3,dqisSig{ggg}];
            numSigVarg3(4,animalInd) = numSigVarg3(4,animalInd) + sum(dqisSig{ggg});
            dK3 = [dK3; dKCoeff_new{ggg}]; dKSig3 = [dKSig3,dkisSig{ggg}];
            numSigVarg3(5,animalInd) = numSigVarg3(5,animalInd) + sum(dkisSig{ggg});
            CKE3 = [CKE3; CKECoeff_new{ggg}]; CKESig3 = [CKESig3,ckeisSig{ggg}];
            numSigVarg3(8,animalInd) = numSigVarg3(8,animalInd) + sum(ckeisSig{ggg});
           
            %end

    end
end

%% plot the fraction of sig grids within each group
% figure;violinplot(numSigVarg1'./size(choiceCoeff{1},1),varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 0.8]);
% ylabel('Fraction of significant grids, group1');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables group1'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group1'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group1'), 'svg');
% 
% figure;violinplot(numSigVarg2'./size(choiceCoeff{1},1),varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 0.8]);
% ylabel('Fraction of significant grids, group2');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables group2'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group2'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group2'), 'svg');
% 
% figure;violinplot(numSigVarg3'./size(choiceCoeff{1},1),varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 0.8]);
% ylabel('Fraction of significant grids, group3');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables group3'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group3'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group3'), 'svg');
% 
% figure;violinplot(numSigVarg4'./size(choiceCoeff{1},1),varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 0.8]);
% ylabel('Fraction of significant grids, group4');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables group4'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group4'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group4'), 'svg');
% 
% % fraction of sig grids (normalized by outcome)
% figure;violinplot((numSigVarg1./numSigVarg1(2,:))',varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 1]);
% ylabel('Fraction of significant grids of outcome, group1');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables of outcome group1'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables of outcome group1'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables of outcome group1'), 'svg');
% 
% figure;violinplot((numSigVarg2./numSigVarg2(2,:))',varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 1]);
% ylabel('Fraction of significant grids of outcome, group2');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables of outcome group2'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables of outcome group2'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables of outcome group2'), 'svg');
% 
% figure;violinplot((numSigVarg3./numSigVarg3(2,:))',varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 1]);
% ylabel('Fraction of significant grids of outcome, group3');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables of outcome group3'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables of outcome group3'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables of outcome group3'), 'svg');
% 
% close all;

%% test for sigficance between different groups
%figure;
%violinplot([numSigVarg1(1,:)./numSigVarg1(2,:);numSigVarg3(1,:)./numSigVarg3(2,:);numSigVarg3(1,:)./numSigVarg3(2,:)]',{'group1','group2','group3'},'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% figure;
% boxplot([numSigVarg1(1,:)./numSigVarg1(2,:);numSigVarg2(1,:)./numSigVarg2(2,:);numSigVarg3(1,:)./numSigVarg3(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(1,:)./numSigVarg1(2,:),numSigVarg2(1,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(1,:)./numSigVarg3(2,:),numSigVarg2(1,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(1,:)./numSigVarg3(2,:),numSigVarg1(1,:)./numSigVarg1(2,:))
% 
% figure;
% boxplot([numSigVarg1(3,:)./numSigVarg1(2,:);numSigVarg2(3,:)./numSigVarg2(2,:);numSigVarg3(3,:)./numSigVarg3(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(3,:)./numSigVarg1(2,:),numSigVarg2(3,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(3,:)./numSigVarg3(2,:),numSigVarg2(3,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(3,:)./numSigVarg3(2,:),numSigVarg1(3,:)./numSigVarg1(2,:))
% 
% figure;
% boxplot([numSigVarg1(4,:)./numSigVarg1(2,:);numSigVarg2(4,:)./numSigVarg2(2,:);numSigVarg3(4,:)./numSigVarg3(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(4,:)./numSigVarg1(2,:),numSigVarg2(4,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(4,:)./numSigVarg3(2,:),numSigVarg2(4,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(4,:)./numSigVarg3(2,:),numSigVarg1(4,:)./numSigVarg1(2,:))
% 
% figure;
% boxplot([numSigVarg1(5,:)./numSigVarg1(2,:);numSigVarg2(5,:)./numSigVarg2(2,:);numSigVarg3(5,:)./numSigVarg3(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(5,:)./numSigVarg1(2,:),numSigVarg2(5,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(5,:)./numSigVarg3(2,:),numSigVarg2(5,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(5,:)./numSigVarg3(2,:),numSigVarg1(5,:)./numSigVarg1(2,:))
% 
% figure;
% boxplot([numSigVarg1(6,:)./numSigVarg1(2,:);numSigVarg2(6,:)./numSigVarg2(2,:);numSigVarg3(6,:)./numSigVarg3(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(6,:)./numSigVarg1(2,:),numSigVarg2(6,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(6,:)./numSigVarg3(2,:),numSigVarg2(6,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(6,:)./numSigVarg3(2,:),numSigVarg1(6,:)./numSigVarg1(2,:))
% 
% figure;
% boxplot([numSigVarg1(7,:)./numSigVarg1(2,:);numSigVarg2(7,:)./numSigVarg2(2,:);numSigVarg3(7,:)./numSigVarg3(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(7,:)./numSigVarg1(2,:),numSigVarg2(7,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(7,:)./numSigVarg3(2,:),numSigVarg2(7,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(7,:)./numSigVarg3(2,:),numSigVarg1(7,:)./numSigVarg1(2,:))
% p = signrank((numSigVarg3(7,:)+numSigVarg2(7,:))./(numSigVarg3(2,:)+numSigVarg2(2,:)),numSigVarg1(7,:)./numSigVarg1(7,:))
% 
% figure;
% boxplot([numSigVarg1(8,:)./numSigVarg1(2,:);numSigVarg2(8,:)./numSigVarg2(2,:);numSigVarg3(8,:)./numSigVarg3(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(8,:)./numSigVarg1(2,:),numSigVarg2(8,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(8,:)./numSigVarg3(2,:),numSigVarg2(8,:)./numSigVarg2(2,:))
% p = signrank(numSigVarg3(8,:)./numSigVarg3(2,:),numSigVarg1(8,:)./numSigVarg1(2,:))
% p = signrank(numSigVarg3(8,:)./numSigVarg3(2,:)+numSigVarg2(8,:)./numSigVarg2(2,:),numSigVarg1(8,:)./numSigVarg1(2,:))
% 
% %% fraction of sig over themselves
% sumAllGroup = numSigVarg1+numSigVarg2+numSigVarg3;
% figure;
% boxplot([[numSigVarg1(1,:);numSigVarg2(1,:);numSigVarg3(1,:)]./sumAllGroup(1,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(1,:)./sumAllGroup(1,:),numSigVarg2(1,:)./sumAllGroup(1,:))
% p = signrank(numSigVarg3(1,:)./sumAllGroup(1,:),numSigVarg2(1,:)./sumAllGroup(1,:))
% p = signrank(numSigVarg3(1,:)./sumAllGroup(1,:),numSigVarg1(1,:)./sumAllGroup(1,:))
% 
% figure;
% boxplot([[numSigVarg1(2,:);numSigVarg2(2,:);numSigVarg3(2,:)]./sumAllGroup(2,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(2,:)./sumAllGroup(2,:),numSigVarg2(2,:)./sumAllGroup(2,:))
% p = signrank(numSigVarg3(2,:)./sumAllGroup(2,:),numSigVarg2(2,:)./sumAllGroup(2,:))
% p = signrank(numSigVarg3(2,:)./sumAllGroup(2,:),numSigVarg1(2,:)./sumAllGroup(2,:))
% 
% figure;
% boxplot([[numSigVarg1(3,:);numSigVarg2(3,:);numSigVarg3(3,:)]./sumAllGroup(3,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(3,:)./sumAllGroup(3,:),numSigVarg2(3,:)./sumAllGroup(3,:))
% p = signrank(numSigVarg3(3,:)./sumAllGroup(3,:),numSigVarg2(3,:)./sumAllGroup(3,:))
% p = signrank(numSigVarg3(3,:)./sumAllGroup(3,:),numSigVarg1(3,:)./sumAllGroup(3,:))
% 
% figure;
% boxplot([[numSigVarg1(4,:);numSigVarg2(4,:);numSigVarg3(4,:)]./sumAllGroup(4,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(4,:)./sumAllGroup(4,:),numSigVarg2(4,:)./sumAllGroup(4,:))
% p = signrank(numSigVarg3(4,:)./sumAllGroup(4,:),numSigVarg2(4,:)./sumAllGroup(4,:))
% p = signrank(numSigVarg3(4,:)./sumAllGroup(4,:),numSigVarg1(4,:)./sumAllGroup(4,:))
% 
% figure;
% boxplot([[numSigVarg1(5,:);numSigVarg2(5,:);numSigVarg3(5,:)]./sumAllGroup(5,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(5,:)./sumAllGroup(5,:),numSigVarg2(5,:)./sumAllGroup(5,:))
% p = signrank(numSigVarg3(5,:)./sumAllGroup(5,:),numSigVarg2(5,:)./sumAllGroup(5,:))
% p = signrank(numSigVarg3(5,:)./sumAllGroup(5,:),numSigVarg1(5,:)./sumAllGroup(5,:))
% 
% figure;
% boxplot([[numSigVarg1(6,:);numSigVarg2(6,:);numSigVarg3(6,:)]./sumAllGroup(6,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(6,:)./sumAllGroup(6,:),numSigVarg2(6,:)./sumAllGroup(6,:))
% p = signrank(numSigVarg3(6,:)./sumAllGroup(6,:),numSigVarg2(6,:)./sumAllGroup(6,:))
% p = signrank(numSigVarg3(6,:)./sumAllGroup(6,:),numSigVarg1(6,:)./sumAllGroup(6,:))
% 
% figure;
% boxplot([[numSigVarg1(7,:);numSigVarg2(7,:);numSigVarg3(7,:)]./sumAllGroup(7,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(7,:)./sumAllGroup(7,:),numSigVarg2(7,:)./sumAllGroup(7,:))
% p = signrank(numSigVarg3(7,:)./sumAllGroup(7,:),numSigVarg2(7,:)./sumAllGroup(7,:))
% p = signrank(numSigVarg3(7,:)./sumAllGroup(7,:),numSigVarg1(7,:)./sumAllGroup(7,:))
% 
% figure;
% boxplot([[numSigVarg1(8,:);numSigVarg2(8,:);numSigVarg3(8,:)]./sumAllGroup(8,:)]','Labels',{'group1','group2','group3'});
% p = signrank(numSigVarg1(8,:)./sumAllGroup(8,:),numSigVarg2(8,:)./sumAllGroup(8,:))
% p = signrank(numSigVarg3(8,:)./sumAllGroup(8,:),numSigVarg2(8,:)./sumAllGroup(8,:))
% p = signrank(numSigVarg3(8,:)./sumAllGroup(8,:),numSigVarg1(8,:)./sumAllGroup(8,:))
% figure;violinplot(numSigVarg4'./size(choiceCoeff{1},1),varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylim([0 0.8]);
% ylabel('Fraction of significant grids, group4');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables group4'));
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group4'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables group4'), 'fig');



% sort group2 and group3
g1.coeff=group1;g1.t = rt;
g1sortOrd = coeff_sort(g1,[0,3]);
g2.coeff =group2; g2.t = rt;
g2sortOrd = coeff_sort(g2,[0,3]);
g3.coeff = group3; g3.t = rt;
g3sortOrd = coeff_sort(g3,[0,3]);

temp1 = group1(g1sortOrd,:); group1 = temp1;
temp2 = group2(g2sortOrd,:); group2 = temp2;
temp3 = group3(g3sortOrd,:); group3 = temp3;

choice1 = choice1(g1sortOrd,:);choice2 = choice2(g2sortOrd,:);choice3 = choice3(g3sortOrd,:);
cSig1 = logical(cSig1(g1sortOrd));cSig2 = logical(cSig2(g2sortOrd)); cSig3 = logical(cSig3(g3sortOrd));
xn1 = xn1(g1sortOrd,:);xn2 = xn2(g2sortOrd,:);xn3 = xn3(g3sortOrd,:);
xSig1 = logical(xSig1(g1sortOrd));xSig2 = logical(xSig2(g2sortOrd)); xSig3 = logical(xSig3(g3sortOrd));
% posRPE1 = posRPE1(g1sortOrd,:);%posRPE2 = posRPE2(g2sortOrd,:);posRPE3 = posRPE3(g3sortOrd,:);
% pRPESig1 = logical(pRPESig1(g1sortOrd));%pRPESig2 = logical(pRPESig2(g2sortOrd)); pRPESig3 = logical(pRPESig3(g3sortOrd));
% negRPE1 = negRPE1(g1sortOrd,:);%negRPE2 = negRPE2(g2sortOrd,:);negRPE3 = negRPE3(g3sortOrd,:);
% nRPESig1 = logical(nRPESig1(g1sortOrd)); %nRPESig2 = logical(nRPESig2(g2sortOrd)); nRPESig3 = logical(nRPESig3(g3sortOrd));
% dQ1 = dQ1(g1sortOrd,:);%dQ2 = dQ2(g2sortOrd,:);dQ3 = dQ3(g3sortOrd,:);
% dQSig1 = dQSig1(g1sortOrd);%dQSig2 = dQSig2(g2sortOrd); dQSig3 = dQSig3(g3sortOrd);
% dK1 = dK1(g1sortOrd,:);%dK2 = dK2(g2sortOrd,:);dK3 = dK3(g3sortOrd,:);
% dKSig1 = dKSig1(g1sortOrd);%dKSig2 = dKSig2(g2sortOrd); dKSig3 = dKSig3(g3sortOrd);
% CKE1 = CKE1(g1sortOrd,:);%CKE2 = CKE2(g2sortOrd,:);CKE3 = CKE3(g3sortOrd,:);
% CKESig1 = CKESig1(g1sortOrd);%CKESig2 = CKESig2(g2sortOrd); CKESig3 = CKESig3(g3sortOrd);

% add one more plot to show the percentage of signficant grids for each
% variables
% venn plot for different groups
% g1Ind = 1:length(cSig1);
% g1cSig = g1Ind(logical(cSig1));g1oSig = g1Ind;g1xSig = g1Ind(logical(xSig1));
% g1dQSig = g1Ind(logical(dQSig1)); g1dKSig = g1Ind(logical(dKSig1)); 
% g1pRPESig = g1Ind(logical(pRPESig1)); g1nRPESig = g1Ind(logical(nRPESig1));
% g1CKESig = g1Ind(logical(CKESig1));
% 
% g2Ind = 1:length(cSig2);
% g2cSig = g2Ind(logical(cSig2));g2oSig = g2Ind;g2xSig = g2Ind(logical(xSig2));
% g2dQSig = g2Ind(logical(dQSig2)); g2dKSig = g2Ind(logical(dKSig2)); 
% g2pRPESig = g2Ind(logical(pRPESig2)); g2nRPESig = g2Ind(logical(nRPESig2));
% g2CKESig = g2Ind(logical(CKESig2));
% 
% g3Ind = 1:length(cSig3);
% g3cSig = g3Ind(logical(cSig3));g3oSig = g3Ind;g3xSig = g3Ind(logical(xSig3));
% g3dQSig = g3Ind(logical(dQSig3)); g3dKSig = g3Ind(logical(dKSig3)); 
% g3pRPESig = g3Ind(logical(pRPESig3)); g3nRPESig = g3Ind(logical(nRPESig3));
% g3CKESig = g3Ind(logical(CKESig3));

% venn plot for all grids

% % c,r,x
% label = {'choice','outcome','interaction'};
% vennPlot(cSigInd,rSigInd,xSigInd,label,savesumfigpath);
% 
% % r-pRPE-nRPE
% label = {'outcome','pRPE','nRPE'};
% vennPlot(rSigInd,pRPESigInd,nRPESigInd,label,savesumfigpath);
% 
% % dQ-pRPE-nRPE
% label = {'dQ','pRPE','nRPE'};
% vennPlot(dQSigInd,pRPESigInd,nRPESigInd,label,savesumfigpath);
% 
% % interaction-pRPE-nRPE
% label = {'interaction','pRPE','nRPE'};
% vennPlot(xSigInd,pRPESigInd,nRPESigInd,label,savesumfigpath);
% 
% % choice-dQ-dK
% label = {'choice','dQ','dK'};
% vennPlot(cSigInd,dQSigInd,dKSigInd,label,savesumfigpath);
% 
% % choice-dK-CKE
% label = {'choice','dK','CKE'};
% vennPlot(cSigInd,dKSigInd,CKESigInd,label,savesumfigpath);

% group1 venn plot

% venn plot for all grids

% % c,r,x
% figure;
% h = vennEulerDiagram({g1cSig;g1oSig;g1xSig}, 'drawProportional', true, 'SetLabels', ["choice"; "outcome"; "interaction"]);
% h.ShowIntersectionCounts = true;
% title('NE group1 c-r-x-dQ-dK significant ROIs');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Ratio of significant grids c-r-x-dq-dk group1'));
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids c-r-x-dq-dk group1'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids c-r-x-dq-dk group1'), 'svg');
% 
% 
% % r-pRPE-nRPE
% figure;
% h = vennEulerDiagram({g1oSig;g1pRPESig;g1nRPESig}, 'drawProportional', true, 'SetLabels', ["outcome"; "pRPE"; "nRPE"]);
% h.ShowIntersectionCounts = true;
% title('NE group1 r-RPE significant ROIs');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Ratio of significant grids r-RPE group1'));
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids r-RPE group1'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids r-RPE group1'), 'svg');
% 
% % dQ-pRPE-nRPE
% figure;
% h = vennEulerDiagram({g1dQSig;g1pRPESig;g1nRPESig}, 'drawProportional', true, 'SetLabels', ["dQ"; "pRPE"; "nRPE"]);
% h.ShowIntersectionCounts = true;
% title('NE group1 dQ-RPE significant ROIs');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Ratio of significant grids dQ-RPE group1'));
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids dQ-RPE group1'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids dQ-RPE group1'), 'svg');
% 
% 
% % interaction-pRPE-nRPE
% figure;
% h = vennEulerDiagram({g1xSig;g1pRPESig;g1nRPESig}, 'drawProportional', true, 'SetLabels', ["interaction"; "pRPE"; "nRPE"]);
% h.ShowIntersectionCounts = true;
% title('NE group1 x-RPE significant ROIs');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Ratio of significant grids x-RPE group1'));
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids x-RPE group1'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids x-RPE group1'), 'svg');
% 
% % choice-dQ-dK
% label = {'g1choice','dQ','dK'};
% vennPlot(g1cSig,g1dQSig,g1dKSig,label,savesumfigpath);
% 
% % choice-dK-CKE
% label = {'g1choice','dK','CKE'};
% vennPlot(g1cSig,g1dKSig,g1CKESig,label,savesumfigpath);
% 
% % group2 venn plot
% figure;
% h = vennEulerDiagram({g2cSig;g2oSig;g2xSig}, 'drawProportional', true, 'SetLabels', ["choice"; "outcome"; "interaction"]);
% h.ShowIntersectionCounts = true;
% title('NE group2 c-r-x significant ROIs');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Ratio of significant grids c-r-x-dq-dk group2'));
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids c-r-x-dq-dk group2'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids c-r-x-dq-dk group2'), 'svg');
% 
% 
% % r-pRPE-nRPE
% figure;
% h = vennEulerDiagram({g2oSig;g2pRPESig;g2nRPESig}, 'drawProportional', true, 'SetLabels', ["outcome"; "pRPE"; "nRPE"]);
% h.ShowIntersectionCounts = true;
% title('NE group2 r-RPE significant ROIs');
% print(gcf,'-dpng',fullfile(savesumfigpath,'Ratio of significant grids r-RPE group2'));
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids r-RPE group2'), 'fig');
% saveas(gcf, fullfile(savesumfigpath,'Ratio of significant grids r-RPE group2'), 'svg');
% 
% % dQ-pRPE-nRPE
% label = {'g2dQ','pRPE','nRPE'};
% vennPlot(g2dQSig,g2pRPESig,g2nRPESig,label,savesumfigpath);
% 
% % interaction-pRPE-nRPE
% label = {'g2interaction','pRPE','nRPE'};
% vennPlot(g2xSig,g2pRPESig,g2nRPESig,label,savesumfigpath);
% 
% % choice-dQ-dK
% label = {'g2choice','dQ','dK'};
% vennPlot(g2cSig,g2dQSig,g2dKSig,label,savesumfigpath);
% 
% % choice-dK-CKE
% label = {'g2choice','dK','CKE'};
% vennPlot(g2cSig,g2dKSig,g2CKESig,label,savesumfigpath);
% 
% % group3 venn plot
% % dQ-pRPE-nRPE
% label = {'g3dQ','pRPE','nRPE'};
% vennPlot(g3dQSig,g3pRPESig,g3nRPESig,label,savesumfigpath);
% 
% % interaction-pRPE-nRPE
% label = {'g3interaction','pRPE','nRPE'};
% vennPlot(g3xSig,g3pRPESig,g3nRPESig,label,savesumfigpath);
% 
% % choice-dQ-dK
% label = {'g3choice','dQ','dK'};
% vennPlot(g3cSig,g3dQSig,g3dKSig,label,savesumfigpath);
% 
% % choice-dK-CKE
% label = {'g3choice','dK','CKE'};
% vennPlot(g3cSig,g3dKSig,g3CKESig,label,savesumfigpath);
% close all;

% outcome
tlabel = 'Outcome';
plot_groupSummary(group1,group2, group3, [], rt, tlabel,savesumfigpath)
tlabel = 'Choice';
plot_groupSummary(choice1,choice2, choice3, [], rt, tlabel,savesumfigpath)
tlabel = 'Interaction';
plot_groupSummary(xn1,xn2, xn3, [], rt, tlabel,savesumfigpath)
tlabel = 'posRPE';
plot_groupSummary(posRPE1,posRPE2, posRPE3, [], rt, tlabel,savesumfigpath)
tlabel = 'negRPE';
plot_groupSummary(negRPE1,negRPE2, negRPE3, [], rt, tlabel,savesumfigpath)
tlabel = 'dQ';
plot_groupSummary(dQ1,dQ2, dQ3, [], rt, tlabel,savesumfigpath)
tlabel = 'dK';
plot_groupSummary(dK1,dK2, dK3, [], rt, tlabel,savesumfigpath)
tlabel = 'CKE';
plot_groupSummary(CKE1,CKE2, CKE3, [], rt, tlabel,savesumfigpath)

tlabel = 'Choice-Sig';
plot_groupSummary(choice1(logical(cSig1),:),choice2(logical(cSig2),:), choice3(logical(cSig3),:), cCoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'Interaction-Sig';
% for interactions, rerun cluster for each grou
plot_groupSummary(xn1(logical(xSig1),:),xn2(logical(xSig2),:), xn3(logical(xSig3),:), xCoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'posRPE-Sig';
plot_groupSummary(posRPE1(logical(pRPESig1),:),posRPE2(logical(pRPESig2),:), posRPE3(logical(pRPESig3),:), posRPECoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'negRPE-Sig';
plot_groupSummary(negRPE1(logical(nRPESig1),:),negRPE2(logical(nRPESig2),:), negRPE3(logical(nRPESig3),:), negRPECoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'dQ-sig';
plot_groupSummary(dQ1(logical(dQSig1),:),dQ2(logical(dQSig2),:), dQ3(logical(dQSig3),:), dQCoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'dK-sig';
plot_groupSummary(dK1(logical(dKSig1),:),dK2(logical(dKSig2),:), dK3(logical(dKSig3),:), dKCoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'CKE-sig';
plot_groupSummary(CKE1(logical(CKESig1),:),CKE2(logical(CKESig2),:), CKE3(logical(CKESig3),:), CKECoeff_nR, rt, tlabel,savesumfigpath)

%% check pos/neg RPE for overlapping space
plot_RPE(posRPECoeff,negRPECoeff,pRPESigInd,nRPESigInd,rt,savesumfigpath)
% check only significant grids (add another group that is not significant
% for outcome but significant for the variable
%% using the average to regroup?
% newGroup = oriGroup;
%
% ref1 = nanmean(group1,1);
% ref2 = nanmean(group2,1);
% ref3 = nanmean(group3,1);
% outInd  =1;
% outlierGroups = cell(0);
%
% for nn = 1:length(rCoeff_new)
%     % ACh, first session as ref group
%
%
%     g1Num = sum(rClust{nn}==1); g2Num = sum(rClust{nn}==2); g3Num = sum(rClust{nn}==3);
%     % calculate the average coefficient with in 1-2.5s for each group
%     g1Mean =nanmean(rCoeff{nn}(rClust{nn}==1,:),1);  g2Mean = nanmean(rCoeff{nn}(rClust{nn}==2,:),1); g3Mean =nanmean(rCoeff{nn}(rClust{nn}==3,:),1);
%
%
%     figure;plot(rt,ref1);hold on;
%     plot(rt,ref2);plot(rt,ref3);
%     plot(rt,g1Mean)
%
%     figure;plot(rt,ref1);hold on;
%     plot(rt,ref2);plot(rt,ref3);
%     plot(rt,g2Mean)
%
%     figure;plot(rt,ref1);hold on;
%     plot(rt,ref2);plot(rt,ref3);
%     plot(rt,g3Mean)
%     corr11 = corrcoef(g1Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr12 = corrcoef(g1Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr13 = corrcoef(g1Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
%     corr21 = corrcoef(g2Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr22 = corrcoef(g2Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr23 = corrcoef(g2Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
%     corr31 = corrcoef(g3Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr32 = corrcoef(g3Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr33 = corrcoef(g3Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
%
%     corr1 = [corr11(1,2),corr12(1,2),corr13(1,2)];
%     corr2 = [corr21(1,2),corr22(1,2),corr23(1,2)];
%     corr3 = [corr31(1,2),corr32(1,2),corr33(1,2)];
%
%     [maxCorr,ind] = max(corr1);
%     if maxCorr>0.5 & g1Num > 25
%         newGroup(nn) = ind;
%     else
%         newGroup(nn) =4; % 4 means outlier
%         outlierGroups{outInd} = rCoeff{nn}(rClust{nn}==1,:);
%         outInd = outInd+1;
%     end
%     [maxCorr,ind] = max(corr2);
%     if maxCorr>0.5 & g2Num > 25
%         newGroup(nn) = ind
%     else
%         newGroup(nn) =4; % 4 means outlier
%         outlierGroups{outInd} = rCoeff{nn}(rClust{nn}==2,:);
%         outInd = outInd+1;
%     end
%     [maxCorr,ind] = max(corr3);
%     if maxCorr>0.5 & g3Num > 25
%         newGroup(nn) = ind;
%
%     else
%         newGroup(nn) =4; % 4 means outlier
%         outlierGroups{outInd} = rCoeff{nn}(rClust{nn}==3,:);
%         outInd = outInd+1;
%     end
% end



% clust.t = rt;
% clust.coeff = clustCoeff;
% tRange = [-3,5]; numClust=6;
% regCoef_cluster(clust,tRange,numClust,tlabel1,'Time from cue(s)',savefluofigpath);


%% get lag
% cnLagCorr = zeros(1,nFiles); cnLagP =  zeros(1,nFiles); cnDis = []; cnCorr = []; cnLag = [];
% rnLagCorr = zeros(1,nFiles); rnLagP = zeros(1,nFiles); rnDis = []; rnCorr = []; rnLag = [];
% xnLagCorr = zeros(1,nFiles); xnLagP = zeros(1,nFiles); xnDis = []; xnCorr = []; xnLag = [];
%
%
% for ii = 1:nFiles
%
%     % load behavior files
%     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
%     %load(fullfile(fn_beh.folder,fn_beh.name));
%
%
%     % load fluorescent files
%     fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
%     if length(fn_fluo) == 1
%
%         savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
%         if ~exist(savefluofigpath,'dir')
%             mkdir(savefluofigpath);
%         end
%
%         cd(savefluofigpath);
%
%
%         %load(fullfile(fn_fluo.folder,fn_fluo.name));
%
%         % make folders to save analysis and plots
%         savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
%         if ~exist(savematpath,'dir')
%             mkdir(savematpath);
%         end
%
%         saveregpath = fullfile(savematpath,'regressionVarTempCorr.mat');  % file path to save the results
%         if exist(saveregpath)
%             load(saveregpath)
%         end
%
%         % get the data
%         % scatter plot of distance from reference grid and lag
%
%         % get distance
%         cnDistance = []; cnlag = []; cnCoeff = [];
%         for tt = 1:size(choicetempData.tempCorrLag,1)
%            for uu = 1:size(choicetempData.tempCorrLag,2)
%                if ~isnan(choicetempData.tempCorrLag(tt,uu))
%                    cnDistance = [cnDistance,sqrt((tt-choicetempData.maxX)^2+(uu-choicetempData.maxY)^2)];
%                    cnlag = [cnlag, choicetempData.tempCorrLag(tt,uu)];
%                    cnCoeff = [cnCoeff, choicetempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%
%         [R,P] = corrcoef(cnDistance,abs(cnlag));
%         if ~isnan(R)
%             cnLagCorr(ii) = R(1,2); cnLagP(ii) = P(1,2); cnDis = [cnDis,cnDistance]; cnCorr = [cnCorr, cnCoeff]; cnLag = [cnLag, cnlag];
%         end
%             %figure;scatter(cnDistance,cnLag);
%
%         rnDistance = []; rnlag = []; rnCoeff = [];
%         for tt = 1:size(outcometempData.tempCorrLag,1)
%            for uu = 1:size(outcometempData.tempCorrLag,2)
%                if ~isnan(outcometempData.tempCorrLag(tt,uu))
%                    rnDistance = [rnDistance,sqrt((tt-outcometempData.maxX)^2+(uu-outcometempData.maxY)^2)];
%                    rnlag = [rnlag, outcometempData.tempCorrLag(tt,uu)];
%                    rnCoeff = [rnCoeff, outcometempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%            [R,P] = corrcoef(rnDistance,abs(rnlag));
%             if ~isnan(R)
%         rnLagCorr(ii) = R(1,2); rnLagP(ii) = P(1,2); rnDis = [rnDis,rnDistance]; rnCorr = [rnCorr, rnCoeff]; rnLag = [rnLag, rnlag];
%             end
%
%         xnDistance = []; xnlag = []; xnCoeff = [];
%         for tt = 1:size(xntempData.tempCorrLag,1)
%            for uu = 1:size(xntempData.tempCorrLag,2)
%                if ~isnan(xntempData.tempCorrLag(tt,uu))
%                    xnDistance = [xnDistance,sqrt((tt-xntempData.maxX)^2+(uu-xntempData.maxY)^2)];
%                    xnlag = [xnlag, xntempData.tempCorrLag(tt,uu)];
%                    xnCoeff = [xnCoeff, xntempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%            [R,P] = corrcoef(xnDistance,abs(xnlag));
%
%             if ~isnan(R)
%         xnLagCorr(ii) = R(1,2); xnLagP(ii) = P(1,2); xnDis = [xnDis,xnDistance]; xnCorr = [xnCorr, xnCoeff]; xnLag = [xnLag, xnlag];
%             end
%     end
% end
%
% % save the results
% savematdir = fullfile(savematsumpath, 'tempCorr');
% if ~exist(savematdir)
%     mkdir(savematdir)
% end
% cnTemp.dis = cnDis; cnTemp.corr = cnCorr; cnTemp.lag = cnLag; cnTemp.corrp = cnLagP; cnTemp.corrCoef = cnLagCorr;
% rnTemp.dis = rnDis; rnTemp.corr = rnCorr; rnTemp.lag = rnLag; rnTemp.corrp = rnLagP; rnTemp.corrCoef = rnLagCorr;
% xnTemp.dis = xnDis; xnTemp.corr = xnCorr; xnTemp.lag = xnLag; xnTemp.corrp = xnLagP; xnTemp.corrCoef = xnLagCorr;
% save(fullfile(savematdir,'tempCorrLag.mat'),'cnTemp','rnTemp','xnTemp');
%
% cnLagCorrSig = cnLagCorr(cnLagP<0.05);
% rnLagCorrSig = rnLagCorr(rnLagP<0.05);
% xnLagCorrSig = xnLagCorr(xnLagP<0.05);
%
%  [p,h,stats] = signrank(cnLagCorrSig,0,'tail','right')
%   [p,h,stats] = signrank(rnLagCorrSig,0,'tail','right')
%    [p,h,stats] = signrank(xnLagCorrSig,0,'tail','right')
%
% figure;scatter(cnDis, abs(cnLag),'.');
% [R,P] = corrcoef(cnDis,abs(cnLag))
%
% figure;scatter(rnDis, abs(rnLag),'.');
% [R,P] = corrcoef(rnDis,abs(rnLag))
%
% figure;scatter(xnDis, abs(xnLag),'.');
% [R,P] = corrcoef(xnDis,abs(xnLag))