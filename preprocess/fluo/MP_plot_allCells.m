function MP_plot_allCells( cells, trialData, savefigpath)
% % plot_allCells %
%PURPOSE:   Plot dF/F for all the cells in one imaging session
%AUTHORS:   MJ Siniscalchi 170105
%           modified by AC Kwan 170515
%
%INPUT ARGUMENTS
%   cells:        Structure generated by calc_dFF().
%   trialData:    Structure generated by flex_getSessionData().
%   blocks:       Structure generated by flex_getBlockData().

%% Display parameters
batchSize=16; %plot in batches of 20 cells

spacing = 2;     % in unit of dF/F
faceAlpha = 0.2; % transparency for rule type colors

nCells = numel(cells.dFF);
t = (cells.t / 60); %Unit: seconds->minutes
startTimes = (trialData.cueTimes ./ 60); %Unit: seconds->minutes

%%
numBatch=ceil(nCells/batchSize);

for l=1:numBatch
    cellidx1=(l-1)*batchSize+1;
    cellidx2=cellidx1+batchSize-1;
    if cellidx2>nCells
        cellidx2=nCells;
    end
    
    %Make color-coded backdrop for each rule block
    figure;
    hold on;
    y1 = spacing*0.5;
    y2 = -[spacing*(batchSize+0.5)]; %Cell idx negated so cell 1 is on top
    %Plot dF/F for all cells
    for i = cellidx1:cellidx2
        if cells.isRedCell{i}
            c = 'r';
        else
            c = 'k';
        end
        plot(t,cells.dFF{i} - spacing*(i-cellidx1+1),'color',c,'LineWidth',1);
    end
    ylabel('dF/F'); xlabel('Time (min)');
    title(['Cells ' int2str(cellidx1) ' to ' int2str(cellidx2) ', out of ' int2str(nCells) '.']);
    xlim([t(1) t(end)]);
    ylim([y2 y1]);
    
    
    print(gcf,'-dpng',fullfile(savefigpath,['dff_allCells-' int2str(l)]));    %png format
    saveas(gcf, fullfile(savefigpath,['dff_allCells-' int2str(l)]), 'fig');
    
end
