function MP_GRAB_tempComp(save_path_mat_NE, save_path_mat_ACh,save_path_fluo_ACh)


varList = {'choice', 'outcome', 'interaction', 'dQ', 'dK', 'posRPE', 'negRPE', 'CKE'};

% get rid of sessions with too few c/x significant grids?

%% compare overlapping area of ACh and NE
sigACh = load(fullfile(save_path_mat_ACh,'sigList.mat'),'SigList');
sigNE = load(fullfile(save_path_mat_NE,'sigList.mat'),'SigList');

nFileACh = max(sigACh.SigList.session);
nFileNE = max(sigNE.SigList.session);

cList_ACh = [];oList_ACh = [];xList_ACh=[];

ocOverlap_ACh = zeros(1,nFileACh);
coOverlap_ACh = zeros(1,nFileACh);
cxOverlap_ACh = zeros(1,nFileACh);
xcOverlap_ACh = zeros(1,nFileACh);
oxOverlap_ACh = zeros(1,nFileACh);
xoOverlap_ACh = zeros(1,nFileACh);

minGrid = 100;
for aa = 1:nFileACh
    c = sigACh.SigList.c(sigACh.SigList.session==aa);
    cList_ACh = [cList_ACh,c];
    o = sigACh.SigList.r(sigACh.SigList.session==aa);
    oList_ACh = [oList_ACh,o];
    x = sigACh.SigList.x(sigACh.SigList.session==aa);
    xList_ACh = [xList_ACh,x];

    
    if sum(c)>minGrid
        ocOverlap_ACh(aa) = sum(c&o)/sum(c);
        coOverlap_ACh(aa) = sum(c&o)/sum(o);
    else
        ocOverlap_ACh(aa)=NaN;
        coOverlap_ACh(aa)=NaN;
    end
      if sum(c)>minGrid & sum(x) > minGrid
            cxOverlap_ACh(aa) = sum(c&x)/sum(x);
            xcOverlap_ACh(aa) = sum(c&x)/sum(c);
    else
        cxOverlap_ACh(aa) = nan;
        xcOverlap_ACh(aa) = nan;
    end
    if sum(x) > minGrid
        oxOverlap_ACh(aa) = sum(o&x)/sum(x);
        xoOverlap_ACh(aa) = sum(o&x)/sum(o);
    else
        oxOverlap_ACh(aa) = NaN;
        xoOverlap_ACh(aa) = NaN;
    end
end
ocOverlap_NE = zeros(1,nFileNE);
coOverlap_NE = zeros(1,nFileNE);
cxOverlap_NE = zeros(1,nFileNE);
xcOverlap_NE = zeros(1,nFileNE);
oxOverlap_NE = zeros(1, nFileNE);
xoOverlap_NE = zeros(1, nFileNE);

cList_NE = [];oList_NE = [];xList_NE=[];


for nn = 1:nFileNE
    c = sigNE.SigList.c(sigNE.SigList.session==nn);
       cList_NE = [cList_NE,c];
    o = sigNE.SigList.r(sigNE.SigList.session==nn);
       oList_NE = [oList_NE,o];
    x = sigNE.SigList.x(sigNE.SigList.session==nn);
        xList_NE = [xList_NE,x];
      
    if sum(c)>minGrid
        ocOverlap_NE(nn) = sum(c&o)/sum(c);
        coOverlap_NE(nn) = sum(c&o)/sum(o);
    else
        coOverlap_NE(nn) = NaN;
        ocOverlap_NE(nn) = NaN;
    end
    if sum(c)>minGrid & sum(x) > minGrid
        cxOverlap_NE(nn) = sum(c&x)/sum(x);
        xcOverlap_NE(nn) = sum(c&x)/sum(c);
    else
        cxOverlap_NE(nn) = nan;
        xcOverlap_NE(nn) = nan;
    end
     if sum(x)>minGrid
      oxOverlap_NE(nn) = sum(o&x)/sum(x);
      xoOverlap_NE(nn) = sum(o&x)/sum(o);
     else
         oxOverlap_NE(nn) = NaN;
         xoOverlap_NE(nn) = NaN;
     end
end

%% chi-square test for independence
%% NE
% choice-outcome
pco = x2testOverlap(cList_NE,oList_NE,length(cList_NE)) 
% choice-interaction
pcx = x2testOverlap(cList_NE,xList_NE,length(cList_NE)) 
% outcome-interaction
pox = x2testOverlap(oList_NE,xList_NE,length(cList_NE)) 

   %% ACh
% choice-outcome
pco = x2testOverlap(cList_ACh,oList_ACh,length(cList_ACh)) 
% choice-interaction
pcx = x2testOverlap(cList_ACh,xList_ACh,length(cList_ACh)) 
% outcome-interaction
pox = x2testOverlap(oList_ACh,xList_ACh,length(cList_ACh)) 


pOC = mediantest(ocOverlap_ACh(~isnan(ocOverlap_ACh)), ocOverlap_NE(~isnan(ocOverlap_NE)))
pCO = mediantest(coOverlap_ACh(~isnan(coOverlap_ACh)), coOverlap_NE(~isnan(coOverlap_NE)))
pOX = mediantest(oxOverlap_ACh(~isnan(oxOverlap_ACh)), oxOverlap_NE(~isnan(oxOverlap_NE)))
pXO = mediantest(xoOverlap_ACh(~isnan(xoOverlap_ACh)), xoOverlap_NE(~isnan(xoOverlap_NE)))
pXC = mediantest(xcOverlap_ACh(~isnan(xcOverlap_ACh)), xcOverlap_NE(~isnan(xcOverlap_NE)))
pCX = mediantest(cxOverlap_ACh(~isnan(cxOverlap_ACh)), cxOverlap_NE(~isnan(cxOverlap_NE)))

pXO = ranksum(oxOverlap_ACh, oxOverlap_NE)
pOC = ranksum(ocOverlap_ACh, ocOverlap_NE,'method','exact')
pCO = ranksum(coOverlap_ACh, coOverlap_NE)
pCX = ranksum(cxOverlap_ACh, cxOverlap_NE)
pXC = ranksum(xcOverlap_ACh, xcOverlap_NE)
pOX = ranksum(oxOverlap_ACh, oxOverlap_NE)
pXO = ranksum(xoOverlap_ACh, xoOverlap_NE)

[h,pOC] = ttest2(ocOverlap_ACh, ocOverlap_NE)
[h,pCO] = ttest2(coOverlap_ACh, coOverlap_NE)
[h,pCX] = ttest2(cxOverlap_ACh, cxOverlap_NE)
[h,pXC] = ttest2(xcOverlap_ACh, xcOverlap_NE)
[h,pXO] = ttest2(xoOverlap_ACh, xoOverlap_NE)
[h,pOX] = ttest2(oxOverlap_ACh, oxOverlap_NE)
varList = {'ACh','NE'};

data.ACh = coOverlap_ACh; data.NE = coOverlap_NE;
figure;violinplot(data,varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
ylabel('Fraction of significant grids');
title('Non-overlap choice over total choice')


Group = [ones(1,length(coOverlap_NE)),2*ones(1,length(coOverlap_ACh))];
figure;
subplot(1,6,1); 
boxplot([ocOverlap_NE,ocOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
ylabel('Fraction of grids');
title('P(r|c)')
ylim([0 1])
subplot(1,6,2); 
boxplot([coOverlap_NE,coOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
ylabel('Fraction of grids');
title('P(c|r)')
ylim([0 1])
subplot(1,6,3)
boxplot([oxOverlap_NE,oxOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(r|x)')
ylim([0 1])
subplot(1,6,4)
boxplot([xoOverlap_NE,xoOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(x|r)')
ylim([0 1])

subplot(1,6,5)
boxplot([cxOverlap_NE,cxOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(c|x)')
ylim([0 1])
subplot(1,6,6)
boxplot([xcOverlap_NE,xcOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(x|c)')
ylim([0 1])
print(gcf,'-dpng',fullfile(save_path_fluo_ACh,'Percentage of overlapping grids with outcome'));
saveas(gcf, fullfile(save_path_fluo_ACh,'Percentage of overlapping grids with outcome'), 'fig');
saveas(gcf, fullfile(save_path_fluo_ACh,'Percentage of overlapping grids with outcome'), 'svg');

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

%% compare the temporal dynamics, outcome
tempNE = fullfile( save_path_mat_NE,'outcome_groupstat.mat');
tempACh = fullfile( save_path_mat_ACh,'outcome_groupstat.mat');
NETemp_o = load(tempNE);
AChTemp_o = load(tempACh);

data1.NE = NETemp_o.outcome_riseT_1;
data1.ACh = AChTemp_o.outcome_riseT_1;
data2.NE = NETemp_o.outcome_riseT_2;
data2.ACh = AChTemp_o.outcome_riseT_2;
data3.NE = NETemp_o.outcome_riseT_3;
data3.ACh = AChTemp_o.outcome_riseT_3;

colors = [63,167,150;255, 189, 53]/255;

figure;
subplot(1,3,1)
violinplot(data1,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
ylabel('Rise time (s)');
title('Group1');
subplot(1,3,2)
violinplot(data2,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group2')
subplot(1,3,3)
violinplot(data3,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group3')
print(gcf,'-dpng',fullfile(savesumfigpath,'rise time outcome var-comp'));
saveas(gcf, fullfile(savesumfigpath,'rise time outcome var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'rise time outcome var-comp'), 'svg');

[h1,p1] = ttest2(NETemp_o.outcome_riseT_1,AChTemp_o.outcome_riseT_1)
[h2,p2] = ttest2(NETemp_o.outcome_riseT_2,AChTemp_o.outcome_riseT_2)
[h3,p3] = ttest2(NETemp_o.outcome_riseT_3,AChTemp_o.outcome_riseT_3)

% report variance

% test for unequal variance
varInd = [ones(1,size(data1.NE,2)),ones(1,size(data1.ACh,2))*2];
p = vartestn([data1.NE,data1.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data2.NE,2)),ones(1,size(data2.ACh,2))*2];
p = vartestn([data2.NE,data2.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data3.NE,2)),ones(1,size(data3.ACh,2))*2];
p = vartestn([data3.NE,data3.ACh]',varInd','TestType','BrownForsythe')

% max time
data1.NE = NETemp_o.outcome_maxT_1;
data1.ACh = AChTemp_o.outcome_maxT_1;
data2.NE = NETemp_o.outcome_maxT_2;
data2.ACh = AChTemp_o.outcome_maxT_2;
data3.NE = NETemp_o.outcome_maxT_3;
data3.ACh = AChTemp_o.outcome_maxT_3;

figure;
subplot(1,3,1)
violinplot(data1,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
ylabel('Rise time (s)');
title('Group1');
subplot(1,3,2)
violinplot(data2,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group2')
subplot(1,3,3)
violinplot(data3,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group3')
print(gcf,'-dpng',fullfile(savesumfigpath,'max time outcome var-comp'));
saveas(gcf, fullfile(savesumfigpath,'max time outcome var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'max time outcome var-comp'), 'svg');

[h1,p1] = ttest2(NETemp_o.outcome_maxT_1,AChTemp_o.outcome_maxT_1)
[h2,p2] = ttest2(NETemp_o.outcome_maxT_2,AChTemp_o.outcome_maxT_2)
[h3,p3] = ttest2(NETemp_o.outcome_maxT_3,AChTemp_o.outcome_maxT_3)

varInd = [ones(1,size(data1.NE,2)),ones(1,size(data1.ACh,2))*2];
p = vartestn([data1.NE,data1.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data2.NE,2)),ones(1,size(data2.ACh,2))*2];
p = vartestn([data2.NE,data2.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data3.NE,2)),ones(1,size(data3.ACh,2))*2];
p = vartestn([data3.NE,data3.ACh]',varInd','TestType','BrownForsythe')
%% compare the temporal dynamics, choice
tempNE = fullfile( save_path_mat_NE,'choice_groupstat.mat');
tempACh = fullfile( save_path_mat_ACh,'choice_groupstat.mat');
NETemp_c = load(tempNE);
AChTemp_c = load(tempACh);

data1.NE = NETemp_c.choice_riseT_1;
data1.ACh = AChTemp_c.choice_riseT_1;
data2.NE = NETemp_c.choice_riseT_2;
data2.ACh = AChTemp_c.choice_riseT_2;
data3.NE = NETemp_c.choice_riseT_3;
data3.ACh = AChTemp_c.choice_riseT_3;

colors = [63,167,150;255, 189, 53]/255;

figure;
subplot(1,3,1)
violinplot(data1,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
ylabel('Rise time (s)');
title('Group1');
subplot(1,3,2)
violinplot(data2,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group2')
ylim([0 5]);
subplot(1,3,3)
violinplot(data3,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group3')
ylim([0 5]);
print(gcf,'-dpng',fullfile(savesumfigpath,'rise time choice var-comp'));
saveas(gcf, fullfile(savesumfigpath,'rise time choice var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'rise time choice var-comp'), 'svg');

[h1,p1] = ttest2(NETemp_c.choice_riseT_1,AChTemp_c.choice_riseT_1)
[h2,p2] = ttest2(NETemp_c.choice_riseT_2,AChTemp_c.choice_riseT_2)
[h3,p3] = ttest2(NETemp_c.choice_riseT_3,AChTemp_c.choice_riseT_3)

varInd = [ones(1,size(data1.NE,2)),ones(1,size(data1.ACh,2))*2];
p = vartestn([data1.NE,data1.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data2.NE,2)),ones(1,size(data2.ACh,2))*2];
p = vartestn([data2.NE,data2.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data3.NE,2)),ones(1,size(data3.ACh,2))*2];
p = vartestn([data3.NE,data3.ACh]',varInd','TestType','BrownForsythe')

% max time
data1.NE = NETemp_c.choice_maxT_1;
data1.ACh = AChTemp_c.choice_maxT_1;
data2.NE = NETemp_c.choice_maxT_2;
data2.ACh = AChTemp_c.choice_maxT_2;
data3.NE = NETemp_c.choice_maxT_3;
data3.ACh = AChTemp_c.choice_maxT_3;

figure;
subplot(1,3,1)
violinplot(data1,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
ylabel('Rise time (s)');
title('Group1');
subplot(1,3,2)
violinplot(data2,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group2')
subplot(1,3,3)
violinplot(data3,varList,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0], 'ShowData' ,false);
title('Group1')
print(gcf,'-dpng',fullfile(savesumfigpath,'max time choice var-comp'));
saveas(gcf, fullfile(savesumfigpath,'max time choice var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'max time choice var-comp'), 'svg');

[h1,p1] = ttest2(NETemp_c.choice_maxT_1,AChTemp_c.choice_maxT_1)
[h2,p2] = ttest2(NETemp_c.choice_maxT_2,AChTemp_c.choice_maxT_2)
[h3,p3] = ttest2(NETemp_c.choice_maxT_3,AChTemp_c.choice_maxT_3)

varInd = [ones(1,size(data1.NE,2)),ones(1,size(data1.ACh,2))*2];
p = vartestn([data1.NE,data1.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data2.NE,2)),ones(1,size(data2.ACh,2))*2];
p = vartestn([data2.NE,data2.ACh]',varInd','TestType','BrownForsythe')
varInd = [ones(1,size(data3.NE,2)),ones(1,size(data3.ACh,2))*2];
p = vartestn([data3.NE,data3.ACh]',varInd','TestType','BrownForsythe')
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



end