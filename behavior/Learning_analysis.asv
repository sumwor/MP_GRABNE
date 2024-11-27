function MP_learning_analysis(dataIndex,save_path)
% % MP_behaviorPerAnimal %
%PURPOSE:   Analyze matching pennies behavior averaged across animals
%AUTHORS:   H Atilgan and AC Kwan 191204
%
%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%   save_path:    path for saving the plots
%
%OUTPUT ARGUMENTS
%

%%
if ~exist(save_path,'dir')
    mkdir(save_path);
end
cd(save_path)
%% go through each animal
% animalFolder = unique(dataIndex.LogFilePath);
% for ii = 1:length(animalFolder)
%     Ind = strfind(animalFolder{ii},filesep);
%     startInd = Ind(end);
%     animalList{ii} = animalFolder{ii}(startInd+1:end);
% end
animalList = unique(dataIndex.Animal);
disp('-----------------------------------------------------------');
disp(['--- Analyzing - summary of ', int2str(numel(animalList)) ' animals']);
disp('-----------------------------------------------------------');

for j = 1:numel(animalList)
    
    %which session belong to this one animal
    currAnimalSessions = ismember(dataIndex.Animal,animalList(j));
    
    %concatenate the sessions for this one animal
    % load the learning result to a cell
    [learning_ses, learning_block] = pennies_learning_summary(dataIndex(currAnimalSessions,:),save_path);
    fieldValues = fieldnames(learning_ses);
    for ff = 1:length(fieldValues)
        if j == 1
            learning_sum_ses.(fieldValues{ff}) = cell(numel(animalList),1);
        end
        learning_sum_ses.(fieldValues{ff}){j}=learning_ses.(fieldValues{ff});
    end

    fieldValues = fieldnames(learning_block);
    for ff = 1:length(fieldValues)
        if j == 1
            learning_sum_block.(fieldValues{ff}) = cell(numel(animalList),1);
        end
        learning_sum_block.(fieldValues{ff}){j}=learning_block.(fieldValues{ff});
    end
end

close all;

% plot the group summary
% align the vectors from different animals to the end of learning
fieldValues = fieldnames(learning_sum_ses);
maxSessions = max(cellfun(@numel, learning_sum_ses.(fieldValues{1})));

for ff = 1:length(fieldValues)
    for aa =1:length(animalList)
        if aa==1
            learning_sum_ses_aligned.(fieldValues{ff}) = nan(maxSessions, length(animalList));
        end
        learning_sum_ses_aligned.(fieldValues{ff})(1:length(learning_sum_ses.(fieldValues{ff}){aa}),aa) = learning_sum_ses.(fieldValues{ff}){aa};
        % end-length(learning_sum_ses.(fieldValues{ff}){aa})+1:end
    end
end

fieldValues = fieldnames(learning_sum_block);
maxBlocks = max(cellfun(@numel, learning_sum_block.(fieldValues{1})));
for ff = 1:length(fieldValues)
    for aa =1:length(animalList)
        if aa==1
            learning_sum_block_aligned.(fieldValues{ff}) = nan(maxBlocks, length(animalList));
        end
        len = length(learning_sum_block.(fieldValues{ff}){aa});
        learning_sum_block_aligned.(fieldValues{ff})(1:len,aa) = flip(learning_sum_block.(fieldValues{ff}){aa});
    end
end

%% generate plot
% bootstrap and make plots
cd(save_path);
figure;

notNaNSessions = sum(~isnan(learning_sum_ses_aligned.(fieldValues{1})),2);
endInd = find(notNaNSessions<3);
%endInd = maxSessions+1;
plotX = 1:endInd(1)-1;
fieldValues = fieldnames(learning_sum_ses);
for ff =1:length(fieldValues)
    plotValue = learning_sum_ses_aligned.(fieldValues{ff})(plotX,:)
    subplot(2,3,ff); hold on;
    plot(plotX,nanmean(plotValue,2),'black')
    hold on;
%errorbar(ITI_interval, meanRT_bin, stdRT_bin);
    ste = nanstd(plotValue,0,2)./sum(~isnan(plotValue),2);
    errorshade(plotX, nanmean(plotValue,2)+ste, nanmean(plotValue,2)-ste,[0.7 0.7 0.7]);
    plot(plotX,nanmean(plotValue,2),'black')
    hold on; 
%     for animal = 1:length(animalList)
%         scatter(plotX,plotValue);
%     end
    ylabel(fieldValues{ff})
    if strcmp(fieldValues{ff},'numTrial')
        ylim([300 800])
    elseif strcmp(fieldValues{ff},'entropy')
        ylim([2 3])
    elseif strcmp(fieldValues{ff},'pmiss')
        ylim([0. 0.3])
    else
        ylim([0.2 0.6])
    end
end
sgtitle(['Learning curve in sessions summary of n=',num2str(length(animalList))])
print(gcf,'-dpng',['Learning curve in sessions summary of n=',num2str(length(animalList))]);    %png format
saveas(gcf, ['Learning curve in sessions summary of n=',num2str(length(animalList))], 'fig');
saveas(gcf, ['Learning curve in sessions summary of n=',num2str(length(animalList))],'svg');

%% block
figure;
fieldValues = fieldnames(learning_sum_block);
notNaNBlocks = sum(~isnan(learning_sum_block_aligned.(fieldValues{1})),2);
endInd = find(notNaNBlocks<3);
%endInd = maxSessions+1;
plotX = 1:endInd(1)-1;
fieldValues = fieldnames(learning_sum_block);
for ff =1:length(fieldValues)
    plotValue = learning_sum_block_aligned.(fieldValues{ff})(plotX,:)
    subplot(2,3,ff); hold on;
    plot(plotX,nanmean(plotValue,2),'black')
    hold on;
%errorbar(ITI_interval, meanRT_bin, stdRT_bin);
    ste = nanstd(plotValue,0,2)./sum(~isnan(plotValue),2);
    errorshade(plotX, nanmean(plotValue,2)+ste, nanmean(plotValue,2)-ste,[0.7 0.7 0.7]);
    plot(plotX,nanmean(plotValue,2),'black')
    hold on; 
%     for animal = 1:length(animalList)
%         scatter(plotX,plotValue);
%     end
    ylabel(fieldValues{ff})
    if strcmp(fieldValues{ff},'numTrial')
        ylim([300 800])
    elseif strcmp(fieldValues{ff},'entropy')
        ylim([2 3])
    elseif strcmp(fieldValues{ff},'pmiss')
        ylim([0. 0.3])
    else
        ylim([0.2 0.6])
    end
end
sgtitle(['Learning curve in sessions summary of n=',num2str(length(animalList))])
print(gcf,'-dpng',['Learning curve in sessions summary of n=',num2str(length(animalList))]);    %png format
saveas(gcf, ['Learning curve in sessions summary of n=',num2str(length(animalList))], 'fig');
saveas(gcf, ['Learning curve in sessions summary of n=',num2str(length(animalList))],'svg');


%plot_logreg(lreg_LR,tlabel);
%print(gcf,'-dpng',fullfile(save_path,'logreg_LR'));    
%saveas(gcf, fullfile(save_path,'logreg_LR'), 'fig');
