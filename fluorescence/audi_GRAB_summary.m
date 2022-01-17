function audi_GRAB_summary(dataIndex)

% summarize all subjects
% Ach right before MP recording: 14, 19,22,25

% load PSTH results, calculate the mean (and STE)
nFiles = size(dataIndex,1);

PSTH_ave = [];
for ii = 1:nFiles
     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    savepsthname = fullfile(fn_beh.folder,'psthMat.mat');
    load(savepsthname)
    psth.t = psth_panel(1).sig{1}.t;
    for jj = 1:length(psth_panel)
        PSTH_ave = [PSTH_ave;psth_panel(jj).sig{1}.signal'];
    end
end

% bootstrap to get the average
psth.coeff= PSTH_ave;

% use bootstrp to get coefficient
psth = getBootstrp(psth, 0, 0.05);
figure;
plot(psth.t,psth.coeff_bootave,'k')
hold on;
gray=[0.7 0.7 0.7];
errorshade(psth.t,psth.bootlow,psth.boothigh,gray);
hold on;
plot(psth.t,psth.coeff_bootave,'k')
set(gca,'box','off');
ax1 = gca; 
ax1.YAxis(1).Visible = 'off'; % remove y-axis
ax1.XAxis.Visible = 'off';
hold on;
plot([2.2 2.7], [0.11 0.11],'k-');
hold on;
plot([2.7 2.7],[0.11 0.12],'k-');
hold on;
plot([0 0],[0.09,0.13],'k--','LineWidth',2)