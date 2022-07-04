function MP_GRAB_tempComp(save_path_mat_NE, save_path_mat_ACh,save_path_fluo_ACh)


varList = {'choice', 'outcome', 'interaction', 'dQ', 'dK', 'posRPE', 'negRPE', 'CKE'};

% get rid of sessions with too few c/x significant grids?

%% compare overlapping area of ACh and NE
sigACh = load(fullfile(save_path_mat_ACh,'sigList.mat'),'SigList');
sigNE = load(fullfile(save_path_mat_NE,'sigList.mat'),'SigList');

nFileACh = max(sigACh.SigList.session);
nFileNE = max(sigNE.SigList.session);

coOverlap_ACh = zeros(1,nFileACh);
cxOverlap_ACh = zeros(1,nFileACh);
oxOverlap_ACh = zeros(1,nFileACh);

minGrid = 50;
for aa = 1:nFileACh
    c = sigACh.SigList.c(sigACh.SigList.session==aa);
    o = sigACh.SigList.r(sigACh.SigList.session==aa);
    x = sigACh.SigList.x(sigACh.SigList.session==aa);
    if sum(c)>minGrid
        coOverlap_ACh(aa) = sum(c&o)/sum(c);
    else
        coOverlap_ACh(aa)=NaN;
    end
      if sum(c)>minGrid & sum(x) > minGrid
     cxOverlap_ACh(aa) = sum(c&x)/sum(c|x);
    else
        cxOverlap_ACh(aa) = nan;
    end
    if sum(x) > minGrid
        oxOverlap_ACh(aa) = sum(o&x)/sum(x);
    else
        oxOverlap_ACh(aa) = NaN;
    end
end
coOverlap_NE = zeros(1,nFileNE);
cxOverlap_NE = zeros(1,nFileNE);
oxOverlap_NE = zeros(1, nFileNE);
for nn = 1:nFileNE
    c = sigNE.SigList.c(sigNE.SigList.session==nn);
    o = sigNE.SigList.r(sigNE.SigList.session==nn);
    x = sigNE.SigList.x(sigNE.SigList.session==nn);
    if sum(c)>minGrid
        coOverlap_NE(nn) = sum(c&o)/sum(c);
    else
        coOverlap_NE(nn) = NaN;
    end
    if sum(c)>minGrid & sum(x) > minGrid
     cxOverlap_NE(nn) = sum(c&x)/sum(c|x);
    else
        cxOverlap_NE(nn) = nan;
    end
     if sum(x)>minGrid
      oxOverlap_NE(nn) = sum(o&x)/sum(x);
     else
         oxOverlap_NE(nn) = NaN;
     end
end

pCO = ranksum(coOverlap_ACh, coOverlap_NE,'tail','left')
pCX = ranksum(cxOverlap_ACh, cxOverlap_NE)
pOX = ranksum(oxOverlap_ACh(~isnan(oxOverlap_ACh)), oxOverlap_NE(~isnan(oxOverlap_NE)),'tail','left')

[h,pCO] = ttest2(coOverlap_ACh, coOverlap_NE,'tail','right')
[h,pCX] = ttest2(cxOverlap_ACh, cxOverlap_NE)
[h,pOX] = ttest2(oxOverlap_ACh, oxOverlap_NE,'tail','right')
varList = {'ACh','NE'};
data.ACh = coOverlap_ACh; data.NE = coOverlap_NE;
figure;violinplot(data,varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
ylabel('Fraction of significant grids');
title('Non-overlap choice over total choice')


Group = [ones(1,length(coOverlap_NE)),2*ones(1,length(coOverlap_ACh))];
figure;
subplot(1,3,1); 
boxplot([coOverlap_NE,coOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
ylabel('Fraction of grids');
title('Choice')
ylim([0 1])
subplot(1,3,2)
boxplot([oxOverlap_NE,oxOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('Interaction')
ylim([0 1])
subplot(1,3,3)
boxplot([cxOverlap_NE,cxOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('Choice-interaction')
ylim([0 1])
print(gcf,'-dpng',fullfile(savesumfigpath,'Percentage of overlapping grids with outcome'));
saveas(gcf, fullfile(savesumfigpath,'Percentage of overlapping grids with outcome'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'Percentage of overlapping grids with outcome'), 'svg');

% data.ACh = cxOverlap_ACh; data.NE = cxOverlap_NE;
% figure;violinplot(data,varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylabel('Fraction of significant grids');
% 
% varList = {'ACh','NE'};
% data.ACh = oxOverlap_ACh; data.NE = oxOverlap_NE;
% figure;violinplot(data,varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
% ylabel('Fraction of significant grids');
% title('Non-overlap interaction over total interaction')

%% ACh
sigList_ACh = load(fullfile(save_path_mat_ACh,'sigNum.mat'));
sigList_NE = load(fullfile(save_path_mat_NE,'sigNum.mat'));
boxGroup = [ones(1,size(sigList_ACh.numSigVar,2)),2*ones(1,size(sigList_NE.numSigVar,2))];

figure;
for pp = 1:length(varList)
subplot(2,4,pp)
boxplot([sigList_ACh.numSigVar(pp,:)./sigList_ACh.numSigTotal,sigList_NE.numSigVar(pp,:)./sigList_NE.numSigTotal],boxGroup, ...
    'Labels',{'ACh','NE'},'PlotStyle','compact','LabelOrientation','horizontal', ...
    'Colors', [255, 189, 53;63,167,150]/255);
title(varList{pp})
set(gca,'box','off')
if pp==1 | pp ==5
    ylabel('Fraction of significant grids');
end
p = ranksum(sigList_ACh.numSigVar(pp,:)./sigList_ACh.numSigTotal,sigList_NE.numSigVar(pp,:)./sigList_NE.numSigTotal);
text(0,1.1,num2str(p))
ylim([0 1]);

end

savesumfigpath = fullfile(save_path_fluo_ACh,'ACh_NE_compare');
if ~exist(savesumfigpath)
    mkdir(savesumfigpath)
end
print(gcf,'-dpng',fullfile(savesumfigpath,'Percentage of significant grids var-comp'));
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids var-comp'), 'svg');


% test for significant difference

% plot boxplot/violin plot with ranksum test. (fraction of significant
% grids)
% savematNE = fullfile(save_path_mat_NE, 'tempCorr','tempCorrLag.mat');
% savematACh = fullfile(save_path_mat_ACh, 'tempCorr','tempCorrLag.mat');
% 
% NETemp = load(savematNE);
% AChTemp = load(savematACh);
% 
% % ,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]
% choice.NE = NETemp.cnTemp.corrCoef;
% choice.ACh = AChTemp.cnTemp.corrCoef;
% figure;violinplot(choice);
% 
% outcome.NE = NETemp.rnTemp.corrCoef;
% outcome.ACh = AChTemp.rnTemp.corrCoef;
% figure;violinplot(outcome);
% 
% interaction.NE = NETemp.xnTemp.corrCoef;
% interaction.ACh = AChTemp.xnTemp.corrCoef;
% figure;violinplot(interaction);
% 
% 
% [p,h]=ranksum(NETemp.cnTemp.corrCoef,AChTemp.cnTemp.corrCoef)
% [p,h]=ranksum(NETemp.rnTemp.corrCoef,AChTemp.rnTemp.corrCoef)
% [p,h]=ranksum(NETemp.xnTemp.corrCoef,AChTemp.xnTemp.corrCoef)