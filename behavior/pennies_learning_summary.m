function [learning_ses, learning_block]=pennies_learning_summary(dataIndex,save_path)



close all;

%% setup path and plotting formats

setup_figprop;  %set up default figure plotting parameters
%set(0, 'DefaultFigureRenderer', 'painters');
%% load data file list
% animalFolder = unique(dataIndex.LogFilePath);
% for ii = 1:length(animalFolder)
%     Ind = strfind(animalFolder{ii},filesep);
%     startInd = Ind(end);
%     animalList{ii} = animalFolder{ii}(startInd+1:end);
% end
animalList = unique(dataIndex.Animal);

%% load the data


% subMask = [];
logInd = 1;
for i = 1:length(animalList)   
    currAnimalSessions = ismember(dataIndex.Animal,animalList(i));
    animalData = dataIndex(currAnimalSessions,:);
    numSessions = length(currAnimalSessions);
    p_reward_ses = zeros(1,numSessions);
    p_switch_ses = zeros(1,numSessions);
    p_miss_ses = zeros(1,numSessions);
    entropy_ses = zeros(1,numSessions);
    sessDate = cell(numSessions,1);
    numTrials_ses = zeros(1,numSessions);
    stats_sub.c=[];  %concatenate choices and outcomes across sessions
    stats_sub.r=[];

    % calculate p(reward) p(switch) p(entropy) for every session

    for jj=1:size(animalData,1)

        date=num2str(animalData.DateNumber(jj));
        yy=date(1:2);
        mm=date(3:4);
        dd=date(5:6);
        sessDate{jj} = [yy mm dd];
        load(fullfile(animalData.BehPath{jj},[yy mm dd],'beh_cut.mat'));
       

        % calculate p_reward, p_switch, entropy
        entropy_ses(jj) = entro;
        p_reward_ses(jj) = nansum(stats.r)/sum(~isnan(stats.r));
        p_miss_ses(jj) = sum(isnan(stats.c(:,1)))/length(stats.r);
        choice_withoutmiss = stats.c(~isnan(stats.c(:,1)),1);
        p_switch_ses(jj) = sum(choice_withoutmiss(1:end-1)~=choice_withoutmiss(2:end))/length(choice_withoutmiss);
        numTrials_ses(jj) = length(stats.r);
        % calculate ITI time for every session
        iti_time = zeros(1, length(trialData.cueTimes)-1);
        
        for tt=1:length(trialData.cueTimes)-1
            iti_time(tt) = trialData.cueTimes(tt+1) - trialData.outcomeTimes(tt);
        end
        iti_trueTime{logInd} = iti_time;
        %iti_num(logInd) = length(trialData.itiTimes);
        trial_num(logInd) = length(trialData.cueTimes);
        
        lick_trType_array{logInd}=lick_trType;
    
        iti_array{logInd}=iti_trType;
        respTime_array{logInd}=respTime_trType;
        trueRespTime{logInd} = trialData.rt;
        choiceBySession{logInd} = stats;
        
        missedTrial{logInd} = trials.miss;

        if exist('nlike_array','var')
            fname = fieldnames(nlike);
            for j=1:numel(fname)  %append for each field
                nlike_array.(fname{j})=[nlike_array.(fname{j}); nlike.(fname{j})];
                bic_array.(fname{j})=[bic_array.(fname{j}); bic.(fname{j})];
            end
            fname = fieldnames(fitpar);
            for j=1:numel(fname)  %append for each field
                fitpar_array.(fname{j})=[fitpar_array.(fname{j}); fitpar.(fname{j})];
            end
        else
            nlike_array = nlike;
            bic_array = bic;
            fitpar_array = fitpar;
        end
        
        nTrial_array(logInd)=sum(stats.c(:,1)==-1)+sum(stats.c(:,1)==1);
        entro_array(logInd)=entro;
        rrate_array(logInd)=sum(stats.r==1)/(sum(stats.r==1)+sum(stats.r==0));
        
        stats_sub.c=[stats_sub.c; stats.c];
        stats_sub.r=[stats_sub.r; stats.r];
        
        logInd = logInd + 1;
    end
    stats_all{i} = stats_sub;
    
    close all;
%     clearvars -except i dirs dataIndex ...
%         lick_trType_array lregRCUC_array lregCRInt_array iti_array respTime_array trueRespTime choiceBySession...
%         nlike_array bic_array fitpar_array iti_trueTime...
%         nTrial_array entro_array rrate_array subMask...
%         stats_all, animalList;
end

%% calculate p/entropy based in 200 trial blocks
block_length = 200;
% calculate from end to start

num_blocks = floor(length(stats_sub.r)/block_length);
entropy_block = zeros(num_blocks,1);
p_reward_block =zeros(num_blocks,1);
p_miss_block = zeros(num_blocks,1);
p_switch_block = zeros(num_blocks,1);
n_block = 1:num_blocks;
nTrials = length(stats_sub.r);
% calculation
for bb =1:num_blocks
    startInd = nTrials - block_length*bb+1; 
    endInd = nTrials-block_length*(bb-1);

    choice = stats_sub.c(startInd:endInd,1);
    reward = stats_sub.r(startInd:endInd);

    p_reward_block(bb) = nansum(reward)/sum(~isnan(reward));
    p_miss_block(bb) = sum(isnan(choice))/length(choice);
    choice_withoutmiss = choice(~isnan(choice),1);
    p_switch_block(bb) = sum(choice_withoutmiss(1:end-1)~=choice_withoutmiss(2:end))/length(choice_withoutmiss);

    % calculate entropy
    choiceBack = 3;
    combos = de2bi([0:2^choiceBack-1],choiceBack);
    combos = 2*(combos - 0.5);  %to make it -1 or 1
    
    % classify animal's choice sequence
    nTrial = size(choice,1);
    cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
    for j = choiceBack:nTrial
        c = choice(j-choiceBack+1:j,1)';
        idx = ismember(combos,c,'rows');   %find if it matches one of the combos
        if sum(idx)==1
            cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
        else
            cumuoccur(:,j) = cumuoccur(:,j-1);
        end
    end
    
    % change the order of the sequence
    p = cumuoccur(:,end)/sum(cumuoccur(:,end));
    entro = 0; % in case of p = 0
    for pp = 1:length(p)
        if p(pp) ~= 0 
            entro = entro-p(pp)*log2(p(pp));
        end
    end

    entropy_block(bb) = entro;
end


%%
savesubFolder = fullfile(save_path,animalList{1});
if ~exist(savesubFolder)
    mkdir(savesubFolder)
end
cd(savesubFolder)
% cd(savebehfigpath);
%save('stats_1.mat', 'choiceBySession','entro_array','respTime_array', 'rrate_array', 'subMask');
% tlabel=strcat('Group summary, n=',int2str(numel(iti_array)), ', subject=',dataIndex.Animal{1});
tlabel=strcat('Group summary, n=',int2str(numel(iti_array)), ', subject=', int2str(length(animalList)));

%% plot the number of missed trials per session from backwards
learning_ses.numTrial = numTrials_ses;
learning_ses.entropy = entropy_ses;
learning_ses.preward = p_reward_ses;
learning_ses.pswitch = p_switch_ses;
learning_ses.pmiss = p_miss_ses;

figure;

plotX = 1:numSessions;
% plot numTrials,p_Miss, entropy, p_reward, p_switch
subplot(2,3,1); hold on;
plot(plotX, numTrials_ses,'black');
ylabel('Number of sessions')

subplot(2,3,2); hold on;
plot(plotX, entropy_ses,'black');
ylabel('Entropy')

subplot(2,3,3); hold on;
plot(plotX, p_reward_ses,'black');
ylabel('P reward')
xlabel(['Sessions'])

subplot(2,3,4); hold on;
plot(plotX, p_miss_ses,'black');
ylabel('P miss')
xlabel(['Sessions'])

subplot(2,3,5); hold on;
plot(plotX, p_switch_ses,'black');
ylabel('P switch');
xlabel('Sessions')
sgtitle(['Learning subject:',animalList{1}])

print(gcf,'-dpng',['Learning subject ',animalList{1}]);    %png format
saveas(gcf, ['Learning subject ',animalList{1}], 'fig');
saveas(gcf, ['Learning subject ',animalList{1}],'svg');

%% plot the number of missed trials per 200-trial block from backwards

learning_block.entropy = entropy_block;
learning_block.preward = p_reward_block;
learning_block.pswitch = p_switch_block;
learning_block.pmiss = p_miss_block;

figure;

plotX = n_block;
% plot numTrials,p_Miss, entropy, p_reward, p_switch

subplot(2,2,1); hold on;
ylim([0 3]);
plot(plotX, flip(entropy_block),'black');
ylabel('Entropy')

subplot(2,2,2); hold on;
ylim([0 1]);
plot(plotX, flip(p_reward_block),'black');
ylabel('P reward')

subplot(2,2,3); hold on;
ylim([0 1]);
plot(plotX, flip(p_miss_block),'black');
ylabel('P miss')
xlabel(['200 trial blocks'])

subplot(2,2,4); hold on;
ylim([0 1]);
plot(plotX, flip(p_switch_block),'black');
ylabel('P switch');
xlabel(['200 trial blocks'])
sgtitle(['Learning subject(block):',animalList{1}])

print(gcf,'-dpng',['Learning subject(block) ',animalList{1}]);    %png format
saveas(gcf, ['Learning subject ',animalList{1}], 'fig');
saveas(gcf, ['Learning subject ',animalList{1}],'svg');


% %% plot the previous ITI and response time
% plot_rtITI(iti_trueTime, trueRespTime)
% print(gcf,'-dpng','rt_ITI_cut');    %png format
% saveas(gcf, 'rt_ITI_cut', 'fig');
% saveas(gcf, 'rt_ITI_cut','svg');
% %%
% plot_lickrate_byTrialType(lick_trType_array);
% 
% plot_val_byTrialType(respTime_array);
% print(gcf,'-dpng','rt_byTrialType_cut');    %png format
% saveas(gcf, 'rt_byTrialType_cut', 'fig');
% saveas(gcf, 'rt_byTrialType_cut','svg');
% 
% plot_val_byTrialType(iti_array);
% print(gcf,'-dpng','iti_byTrialType_cut');   %png format
% saveas(gcf, 'iti_byTrialType_cut', 'fig');
% saveas(gcf, 'iti_byTrialType_cut','svg');
% 
% %% plot iti_array in tertials
% iti = [];
% for ii = 1:length(iti_trueTime)
%     iti = [iti,iti_trueTime{ii}];
% end
% shortEdge = 5.26;
% mediumEdge = 6.50;
% 
% figure;
% edges1 = 4:0.1:round(shortEdge*10)/10;
% histogram(iti(iti<shortEdge),edges1,'EdgeColor','none','FaceColor',[0 0.4470 0.7410],'FaceAlpha',1);
% hold on;
% edges2 = round(shortEdge*10)/10:0.1:round(mediumEdge*10)/10;
% histogram(iti(iti<mediumEdge & iti>=shortEdge),edges2,'EdgeColor','none','FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',1);
% hold on;
% edges3 = round(mediumEdge*10)/10:0.1:max(iti);
% histogram(iti(iti>mediumEdge),edges3,'EdgeColor','none','FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',1);
% set(gca,'box','off');
% print(gcf,'-dpng','iti_byLength');   %png format
% saveas(gcf, 'iti_byLength', 'fig');
% saveas(gcf, 'iti_byLength','svg');
% 
% %%
% %disp(['Total number of trials: ' int2str(sum(~isnan(stats_all.c(:,1))))]);
% disp(['Model fit, using concatenated choice behaviors:']);
% 
% model{1}.name = 'WSLS';
%     model{1}.fun = 'funWSLS';
%     model{1}.initpar=0.5; % initial [prob_WSLS]
%     model{1}.lb=0;
%     model{1}.ub=1;
% 
%     
%     model{2}.name = 'Q_RPE';
%     model{2}.fun = 'funQ_RPE';
%     model{2}.initpar=[0.5 10]; % initial [alpha beta]
%     model{2}.lb=[0 0];
%     model{2}.ub=[1 inf];
% 
%     
%     model{3}.name = 'DQ_RPE';           % text label to refer to the model
%     model{3}.fun = 'funDQ_RPE';     % the corresponding .m code for the model
%     model{3}.initpar=[0.5 5 0.2];   % initial [alpha_reward beta alpha_noreward]
%     model{3}.lb=[0 0 0];            % upper bound of parameters
%     model{3}.ub=[1 inf 1];          % lower bound of parameters
% 
% 
%     model{4}.name = 'FQ_RPE';      % text label to refer to the model
%     model{4}.fun = 'funFQ_RPE';    % the corresponding .m code for the model
%     model{4}.initpar=[0.5 5];      % initial [alpha_reward beta]
%     model{4}.lb=[0 0];             % upper bound of parameters
%     model{4}.ub=[1 inf];           % lower bound of parameters
%     
%     model{5}.name = 'FQ_RPE_CK'; % with a choice autocorrelation term
%     model{5}.fun = 'funFQ_RPE_CK';
%     model{5}.initpar = [0.1 1 0 0]; % initial [alpha beta tau phi]
%     model{5}.lb = [0 0 0 0];
%     model{5}.ub = [1 inf 1 inf];
%    
%     
%     
%  for hh = 1:length(animalList)
%     stats_fit.c = stats_all{hh}.c(:,1);
%     stats_fit.r = stats_all{hh}.r;
%     fitpar = struct;
%     bic = struct;
%     nlike = struct;
% %     fitpar = cell(0);
% %     bic = cell(0);
%     for kk=1:5
%         if isfield(model{kk},'lb')
%             [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar,model{kk}.lb,model{kk}.ub);
%         else
%             [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar);
%         end
%     end
%     Fitpar{hh} = fitpar;
%     Bic{hh} = bic;
%     Nlike{hh} = nlike;
%  end
% 
% %[~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats_all,1,2);
% %[~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats_all,1,5);
% %[~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats_all,1,10);
% 
% 
% %% simulate to get latent action value, compare with animal's choice
% player1.label='algo_FQ_RPE';   % change to CA later
% player1.params.a=Fitpar{1}.FQ_RPE(1);
% player1.params.b=Fitpar{1}.FQ_RPE(2);
% 
% stats_sim=predictAgent(player1,stats_all{1});
% 
% x = 1;  %player 1
% n_plot = 500;   %plot first 500 trials
% plot_session_qparam(stats_sim,x,n_plot);
% for jj = 1:length(animalList)
%     player1.params.a=Fitpar{jj}.FQ_RPE(1);
%     player1.params.b=Fitpar{jj}.FQ_RPE(2);
% 
%     stats_simAll{jj} = predictAgent(player1,stats_all{jj});
% end
% plot_session_qhist(stats_simAll,x)
% 
% 
% player1.label='algo_FQ_RPE_CK';   % change to CA later
% player1.params.a=Fitpar{1}.FQ_RPE_CK(1);
% player1.params.b=Fitpar{1}.FQ_RPE_CK(2);
% player1.params.ac = Fitpar{1}.FQ_RPE_CK(3);
% player1.params.bc = Fitpar{1}.FQ_RPE_CK(4);
% stats_sim=predictAgent(player1,stats_all{1});
% 
% x = 1;  %player 1
% n_plot = 500;   %plot first 500 trials
% plot_session_qparam(stats_sim,x,n_plot);
% for jj = 1:length(animalList)
%     player1.params.a=Fitpar{jj}.FQ_RPE_CK(1);
%     player1.params.b=Fitpar{jj}.FQ_RPE_CK(2);
%     player1.params.ac = Fitpar{jj}.FQ_RPE_CK(3);
%     player1.params.bc = Fitpar{jj}.FQ_RPE_CK(4);
%     stats_simAll{jj} = predictAgent(player1,stats_all{jj});
% end
% plot_session_qhist(stats_simAll,x)
% 
% 
% %%
% disp('---Mean normalized likelihood for fits per session');
% fname = fieldnames(nlike_array);
% for j=1:numel(fname)  %append for each field
%     disp([fname{j} ' - ' num2str(nanmean(nlike_array.(fname{j})))]);
% end
% 
% disp('---Mean BIC for fits per session');
% fname = fieldnames(bic_array);
% for j=1:numel(fname)  %append for each field
%     disp([fname{j} ' - ' num2str(nanmean(bic_array.(fname{j})))]);
% end
% 
% %%
% figure;
% subplot(2,3,1); hold on;
% %plot(rand(1,numel(nTrial_array)),nTrial_array,'k^','MarkerSize',15);
% %boxplot(nTrial_array,'Colors','k','Notch','off','Labels',[char(956),'=',num2str(round(mean(nTrial_array)))]);
% boxplot(nTrial_array,'Colors','k','Notch','off');
% %boxplot(nTrial_array,'Colors','k','Notch','off');
% ylim([0 1000]); 
% ylabel(['Trials performed']);
% set(gca,'box','off') 
% 
% subplot(2,3,2); hold on;
% %plot(rand(1,numel(entro_array)),entro_array,'k^','MarkerSize',15);
% %boxplot(entro_array,'Colors','k','Symbol','k+','Notch','off','Labels',{[char(956),'=',num2str((mean(entro_array)))]});
% boxplot(entro_array,'Colors','k','Symbol','k+','Notch','off');
% 
% plot([-1 2],[3 3],'k--','LineWidth',2);
% ylim([2 3.1]); 
% ylabel('Entropy (bits)');
% set(gca,'box','off') 
% 
% subplot(2,3,3); hold on;
% %plot(rand(1,numel(rrate_array)),100*rrate_array,'k^','MarkerSize',15);
% %boxplot(rrate_array*100,'Colors','k','Symbol','k+','Notch','off','Labels',{[char(956),'=',num2str((mean(rrate_array))*100),'%']});
% boxplot(rrate_array*100,'Colors','k','Symbol','k+','Notch','off');
% plot([-1 2],[50 50],'k--','LineWidth',2);
% ylim([20 65]);
% ylabel('Reward rate (%)');
% set(gca,'box','off') 
% print(gcf,'-dpng',['summary' int2str(x)]);    %png format
% saveas(gcf,['summary' int2str(x)], 'fig');
% saveas(gcf, ['summary' int2str(x)],'svg');
% 
% figure;
% subplot(1,3,1); hold on;
% plot([-0.1 1.1],[0 0],'k--','LineWidth',2);
% plot([0 0],[-0.4 4],'k--','LineWidth',2);
% plot(fitpar_array.FQ_RPE(:,1),fitpar_array.FQ_RPE(:,2),'k.','MarkerSize',30);
% xlabel('\alpha, Learning rate');
% ylabel('\beta, Inverse temperature');
% xlim([-0.1 1.1]);
% ylim([-0.4 4]);
% 
% subplot(1,3,2); hold on;
% plot([-0.1 1.1],[0 0],'k--','LineWidth',2);
% plot([0 0],[-0.4 4],'k--','LineWidth',2);
% plot(fitpar_array.FQ_RPE_CK(:,1),fitpar_array.FQ_RPE_CK(:,2),'k.','MarkerSize',30);
% xlabel('\alpha, Learning rate (CK)');
% ylabel('\beta, Inverse temperature (CK)');
% xlim([-0.1 1.1]);
% ylim([-0.4 4]);
% 
% subplot(1,3,3); hold on;
% plot([-0.1 1.1],[0 0],'k--','LineWidth',2);
% plot([0 0],[-0.4 4],'k--','LineWidth',2);
% plot(fitpar_array.FQ_RPE_CK(:,3),fitpar_array.FQ_RPE_CK(:,4),'k.','MarkerSize',30);
% xlabel('\alpha_c, Learning rate (CK)');
% ylabel('\beta_c, Inverse temperature (CK)');
% xlim([-0.1 1.1]);
% ylim([-0.4 4]);
% print(gcf,'-dpng',['alpha-beta_cut' int2str(x)]);    %png format
% saveas(gcf,['alpha-beta_cut' int2str(x)], 'fig');
% saveas(gcf, ['alpha-beta_cut' int2str(x)],'svg');
% % 
% % alpha = fitpar_array.Q_RPE(:,1); beta = fitpar_array.Q_RPE(:,2);
% % save('a_b_saline.mat','alpha','beta');
% 
% %% plot the alphas and betas and relative value of beta/beta_k
% beta = zeros(1,length(animalList));
% beta_K = zeros(1,length(animalList));
% for ii = 1:length(animalList)
%     beta(ii) = Fitpar{ii}.FQ_RPE_CK(2);
%     beta_K(ii) = Fitpar{ii}.FQ_RPE_CK(4);
% end
% rel_beta = beta_K./beta;
% 
% figure;
% h = boxplot(rel_beta,'Colors','k','Symbol','k+','Notch','off');
% ylabel('\beta_K/\beta','FontSize',50)
% ylim([0 inf]);
% %set(gca,'YTick',[0:2:inf]);
% set(h,{'linew'},{4})
% set(gca,'box','off');
% set(gca,'linewidth',4);
% print(gcf,'-dpng',['betaK-beta' int2str(x)]);    %png format
% saveas(gcf,['betaK-beta' int2str(x)], 'fig');
% saveas(gcf, ['betaK-beta' int2str(x)],'svg');
% 
% alpha = zeros(1,length(animalList));
% alpha_K = zeros(1,length(animalList));
% for ii = 1:length(animalList)
%     alpha(ii) = Fitpar{ii}.FQ_RPE_CK(1);
%     alpha_K(ii) = Fitpar{ii}.FQ_RPE_CK(3);
% end
% rel_alpha = alpha_K./alpha;
% 
% % alpha/alpha_K
% figure;
% h = boxplot(rel_alpha,'Colors','k','Symbol','k+','Notch','off');
% ylabel('\alpha_K/\alpha','FontSize',50)
% %ax=gca;
% %ax.YAxisLocation = 'origin';
% ylim([0 1.6]);
% %set(gca,'YTick',[0:0.2:1.2]);
% set(h,{'linew'},{4})
% set(gca,'box','off');
% set(gca,'linewidth',4);
% print(gcf,'-dpng',['alphaK-alpha' int2str(x)]);    %png format
% saveas(gcf,['alphaK-alpha' int2str(x)], 'fig');
% saveas(gcf, ['alphaK-alpha' int2str(x)],'svg');
% 
% % beta_K + beta
% sum_beta = beta_K + beta;
% 
% figure;
% h = boxplot(sum_beta,'Colors','k','Symbol','k+','Notch','off');
% ylabel('\beta_K+\beta','FontSize',50)
% ylim([0 6]);
% %set(gca,'YTick',[0:1:4]);
% set(h,{'linew'},{4})
% set(gca,'box','off');
% set(gca,'linewidth',4);
% print(gcf,'-dpng',['beta_K+beta' int2str(x)]);    %png format
% saveas(gcf,['beta_K+beta' int2str(x)], 'fig');
% saveas(gcf, ['beta_K+beta' int2str(x)],'svg');
% 
% % beta_K / (beta_K + beta)
% frac_beta = beta_K./(sum_beta);
% figure;
% h = boxplot(frac_beta,'Colors','k','Symbol','k+','Notch','off');
% ylabel('\beta_K/(\beta_K+\beta)','FontSize',50)
% ylim([0 1]);
% set(gca,'YTick',[0:0.2:1]);
% set(h,{'linew'},{4})
% set(gca,'box','off');
% set(gca,'linewidth',4);
% print(gcf,'-dpng',['frac_beta' int2str(x)]);    %png format
% saveas(gcf,['frac_beta' int2str(x)], 'fig');
% saveas(gcf, ['frac_beta' int2str(x)],'svg');
