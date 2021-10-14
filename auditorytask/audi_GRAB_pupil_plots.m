function audi_GRAB_pupil_plots(dataIndex)

% load behavior and dff data
% make some simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    % load dFF files
    load(fullfile(fn_beh.folder,'dff.mat'));
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii},['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    load(fullfile(fn_pup.folder,fn_pup.name));
    
    savefigpath = fullfile(fn_beh.folder,'figs-fluo');
    if ~exist(savefigpath,'dir')
        mkdir(savefigpath);
    end
    cd(savefigpath);

        %% Plot dF/F of all the cells

        %% get cue-aligned dF/F for each cell
        % Fig. 3d in paper came from 140605 data set, cell 8 10 37 74
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-1:0.05:3];
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = 0.95;  %confidence interval
        params.minNumTrial = 50;
        for j=1:numel(cells.dFF)
            psth_output=[];
            
            psth_panel(j).sig{1} = get_psth( cells.dFF{j}, cells.t, params.trigTime(2:end), 'df/f', params );
           
            tlabel = ['Cell #',num2str(j)];
            %plot_psth(psth_panel,tlabel,params.xtitle);
            %print(gcf,'-dpng',fullfile(savefigpath,['cell' int2str(j)]));
            %close;
        end
        params.window = [-1:0.1:3];
        psthpupil_panel(1).sig{1} = get_psth(pupil.dia, pupil.t, params.trigTime(2:end-1), 'Pupil change', params );
          
        gray=[0.7 0.7 0.7];
        avePsth = zeros(size(psth_panel(1).sig{1}.signal,1),1);
        figure;
        for tt=1:length(cells.dFF)
            hold on; 
            t=psth_panel(tt).sig{1}.t;

             plot(t,psth_panel(tt).sig{1}.signal,'Color',[0.8500 0.3250 0.0980 0.2],'LineWidth',0.5);
             avePsth = avePsth + psth_panel(tt).sig{1}.signal/length(cells.dFF);
             %errorshade(t,psth_panel(tt).sig{1}.bootlow,psth_panel(tt).sig{1}.boothigh,gray);
        end
        ax1 = gca; 
        ax1.YAxis(1).Visible = 'off'; % remove y-axis
        ax1.XAxis.Visible = 'off';
        hold on; plot(t, avePsth,'-','Color',[0.8500 0.3250 0.0980]);
        ylim([0.05 0.25]);
        % get average psth
        yyaxis right
        hold on; 
        
       errorshade(psthpupil_panel(1).sig{1}.t,psthpupil_panel(1).sig{1}.bootlow,psthpupil_panel(1).sig{1}.boothigh,[0 0.4470 0.7410], 0.2);
        % 0.95 CI is too large
       hold on; plot(psthpupil_panel(1).sig{1}.t,psthpupil_panel(1).sig{1}.signal,'-','Color',[0 0.4470 0.7410]);
       hold on; 
       ylim([-0.3 0.6]);
       ax1 = gca; 
       ax1.YAxis(2).Visible = 'off'; % remove y-axis
       ax1.XAxis.Visible = 'off'; % remove y-axis
        set(gca,'box','off');
        
         hold on;
             plot([2.5 2.7], [0.18 0.18],'k-');
             hold on;
             plot([2.7 2.7],[0.18 0.2],'k-');
             hold on;
             plot([0 0],[-0.5 0.5],'k--');
                print(gcf,'-dpng',['tone-fluo_pupil_ave' ]);
             saveas(gcf, 'tone-fluo_pupil_ave', 'fig');
             saveas(gcf, 'tone-fluo_ave', 'svg');
end