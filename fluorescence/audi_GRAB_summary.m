function audi_GRAB_summary(dataIndex,savefigpath)

% summarize all subjects
% Ach right before MP recording: 14, 19,22,25
% NE tone response
% load PSTH results, calculate the mean (and STE)

Ach = [19,27:35];
NE = [20:26];
nFiles = size(dataIndex,1);

PSTH_ave_Ach = [];
PSTH_ave_AchPup = [];
PSTH_ave_NE = [];
PSTH_ave_NEPup= [];
for ii = 1:length(Ach)
     fn_beh = dir(fullfile(dataIndex.BehPath{Ach(ii)},'*beh.mat'));
      date = num2str(dataIndex.DateNumber(Ach(ii)));
    fn_pup = dir(fullfile(dataIndex.BehPath{Ach(ii)},['*',date(1:6),'*_pup.mat']));
    savepsthname = fullfile(fn_beh.folder,'psthMat.mat');
      savepupname = fullfile(fn_pup.folder,'psthPup.mat');
    load(savepsthname);
    load(savepupname);
    psth_Ach.t = psth_panel(1).sig{1}.t;
    psth_AchPup.t = psthpupil_panel(1).sig{1}.t;
    for jj = 1:length(psth_panel)
        PSTH_ave_Ach = [PSTH_ave_Ach;psth_panel(jj).sig{1}.signal'];
    end
    PSTH_ave_AchPup = [PSTH_ave_AchPup;psthpupil_panel(1).sig{1}.signal'];
end

for ii = 1:length(NE)
     fn_beh = dir(fullfile(dataIndex.BehPath{NE(ii)},'*beh.mat'));
      date = num2str(dataIndex.DateNumber(NE(ii)));
    fn_pup = dir(fullfile(dataIndex.BehPath{NE(ii)},['*',date(1:6),'*_pup.mat']));
    savepsthname = fullfile(fn_beh.folder,'psthMat.mat');
      savepupname = fullfile(fn_pup.folder,'psthPup.mat');
    load(savepsthname)
      load(savepupname);
    psth_NE.t = psth_panel(1).sig{1}.t;
     psth_NEPup.t = psthpupil_panel(1).sig{1}.t;
    for jj = 1:length(psth_panel)
        PSTH_ave_NE = [PSTH_ave_NE;psth_panel(jj).sig{1}.signal'];
    end
    PSTH_ave_NEPup = [PSTH_ave_NEPup;psthpupil_panel(1).sig{1}.signal'];
end

% bootstrap to get the average
psth_Ach.coeff= PSTH_ave_Ach;
psth_NE.coeff= PSTH_ave_NE;
psth_AchPup.coeff = PSTH_ave_AchPup;
psth_NEPup.coeff = PSTH_ave_NEPup;
% use bootstrp to get coefficient
psth_Ach = getBootstrp(psth_Ach, 0, 0.05);
psth_NE = getBootstrp(psth_NE, 0, 0.05);
psth_AchPup = getBootstrp(psth_AchPup, 0, 0.05);
psth_NEPup = getBootstrp(psth_NEPup, 0, 0.05);
% normalize baseline
baseAch = nanmean(psth_Ach.coeff_bootave(psth_Ach.t<0));
baseNE = nanmean(psth_NE.coeff_bootave(psth_NE.t<0));
baseAchPup = nanmean(psth_AchPup.coeff_bootave(psth_AchPup.t<0));
baseNEPup = nanmean(psth_NEPup.coeff_bootave(psth_NEPup.t<0));

figure;
plot(psth_Ach.t,psth_Ach.coeff_bootave-baseAch,'Color',[255 189 53]/255)
hold on;
plot(psth_NE.t,psth_NE.coeff_bootave-baseNE,'Color',[63,167,150]/255)
hold on;
errorshade(psth_Ach.t,psth_Ach.bootlow-baseAch,psth_Ach.boothigh-baseAch,[255 189 53]/255,0.5);
hold on;
errorshade(psth_NE.t,psth_NE.bootlow-baseNE,psth_NE.boothigh-baseNE,[63,167,150]/255,0.5);
hold on;
plot(psth_NE.t,psth_NE.coeff_bootave-baseNE,'Color',[63,167,150]/255)
hold on;
plot(psth_Ach.t,psth_Ach.coeff_bootave-baseAch,'Color',[255,189,53]/255)

% plot pupil
% hold on;plot(psth_AchPup.t, psth_AchPup.coeff_bootave-baseAchPup);
% hold on;plot(psth_NEPup.t, psth_NEPup.coeff_bootave-baseNEPup);
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('df/f');
legend('Ach','NE');
ax1 = gca; 


title('GRAB tone response');
print(gcf,'-dpng',fullfile(savefigpath,'GRAB-tone'));    %png format
saveas(gcf, fullfile(savefigpath,'GRAB-tone'), 'fig');
saveas(gcf, fullfile(savefigpath,'GRAB-tone'),'svg');
