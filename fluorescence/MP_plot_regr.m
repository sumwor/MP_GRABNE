function MP_plot_regr(input, input_ctrl, pvalThresh,tlabel,xtitle,control)
% % plot_regr %
%PURPOSE:   Plot results from multiple linear regression
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   input:        Structure generated by linear_regr().
%   input_ctrl:   control values of significant sessions, if [] indicate no
%                   contrl, then use binomial test
%   pvalThresh:   Threshold value to deem whether p-values are significant
%   tlabel:       Text to put as title of the plot.
%   xtitle:       Text to put as the label for x-axis.

%% setup

% check whether a control regression exist
if ~exist('control','var')
    % third parameter does not exist, so default it to something
    exist_control = 0;
else
    exist_control = 1;
end

t=input{1}.regr_time;
dt=nanmean(diff(t));

nCells=numel(input);
for j=1:nCells
    pval(:,:,j)=input{j}.pval;
end

nPredictor=input{1}.numPredictor;

nback=input{1}.nback;

if (input{1}.interaction == true)
    if nPredictor == 2
        panelv = nPredictor + 1;    %plot extra row for the interaction terms
        nInteraction = 1;
    else
        panelv = nPredictor + 2;    %plot extra row for the interaction terms
        nInteraction = 2;
    end
else
    panelv = nPredictor;
    nInteraction = 0;
end

%% plot results

figure;

for l=1:nPredictor
    for k=1:1+nback
        currPredictor=1+k+(l-1)*(1+nback); %first +1 because first term is bias
        
        if nPredictor == 7  % C(n+1) regression
            subplot(2,4,currPredictor-1); hold on;
             [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 5
            subplot(2,3,currPredictor-1); hold on;
             [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 2 & nback == 0  % Pos/neg RPE
            subplot(1,2,currPredictor-1); hold on;
             [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 14 % for future choice and reward
            subplot(4,4, currPredictor-1); hold on;
             [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 10  % for dQ and chosenQ plot
            subplot(4,3,currPredictor-1); hold on;
             [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 6  % for RPE/updatedQ
            subplot(2,3,currPredictor-1); hold on;
        elseif nPredictor == 12  % for RL with choice autocorrelation
            subplot(4,3,currPredictor-1); hold on;
             [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 9  % for RL with choice autocorrelation
            subplot(3,3,currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        else
            subplot(panelv,1+nback,currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback)
        end
        %patch([t(1) t(end) t(end) t(1)],[0 0 100*pvalThresh 100*pvalThresh],[0.7 0.7 0.7],'EdgeColor','none');
        percentage = 100*sum(pval(:,currPredictor,:)<pvalThresh,3)/nCells;
        %plot(t,percentage,'k.-','MarkerSize',15);
        %if isempty(input_ctrl)
        plot(t,percentage,'k.-');
%         else
%             fieldName = ['trigEvent',num2str(currPredictor-1)];
%             for ll=1:numel(percentage)
%                 if input_ctrl.(fieldName).boothigh(ll) < percentage(ll)
%                     if ll==1
%                         plot(t(ll)+dt*[-0.5 0.5],[(percentage(ll)+percentage(ll))/2 (percentage(ll)+percentage(ll+1))/2],'-','Color',[0 0 0],'LineWidth',5);
%                     elseif ll==length(percentage)
%                         plot(t(ll)+dt*[-0.5 0.5],[(percentage(ll)+percentage(ll))/2 (percentage(ll)+percentage(ll))/2],'g-','LineWidth',5);
%                     else
%                         plot(t(ll)+dt*[-0.5 0],[(percentage(ll)+percentage(ll-1))/2 percentage(ll)],'-','Color',[0 0 0],'LineWidth',7);
%                         plot(t(ll)+dt*[0 0.5], [percentage(ll) (percentage(ll)+percentage(ll+1))/2],'-','Color',[0 0 0],'LineWidth',7);
%                     end
%                 else % not significant
%                     if ll==1
%                         plot(t(ll)+dt*[-0.5 0.5],[(percentage(ll)+percentage(ll))/2 (percentage(ll)+percentage(ll+1))/2],'k-','LineWidth',5);
%                     elseif ll==length(percentage)
%                         plot(t(ll)+dt*[-0.5 0.5],[(percentage(ll)+percentage(ll))/2 (percentage(ll)+percentage(ll))/2],'k-','LineWidth',5);
%                     else
%                         plot(t(ll)+dt*[-0.5 0],[(percentage(ll)+percentage(ll-1))/2 percentage(ll)],'k-','LineWidth',7);
%                         plot(t(ll)+dt*[0 0.5], [percentage(ll) (percentage(ll)+percentage(ll+1))/2],'k-','LineWidth',7);
%                     end
%                     
%                 end
%             end
%         end
        plot([0 0],[0 100],'k','LineWidth',1);
        xlim([floor(t(1)) ceil(t(end))]);
        xticks([floor(t(1)):1:ceil(t(end))]);
        ylim([0 104]);
        title(tlabel{currPredictor-1});
        if if_xlabel ~= 1
            set(gca,'xticklabel',[])
        end
        if if_xlabel == 1
            xlabel(xtitle);
        end
        if if_ylabel == 1
            ylabel_txt = {'Fraction of','sessions (%)'};
            ylabel(ylabel_txt);
        end

        %identifying significant points via binomial test
        if isempty(input_ctrl) % if no input control regression, use binomial test
            sig=[];
            for ll=1:numel(t)
                [p]=myBinomTest(sum(pval(ll,currPredictor,:)<pvalThresh,3),nCells,pvalThresh);
                sig(ll)=p;
            end
            %ifSig = sig<pvalThresh/10;
            %start1 = strfind([0, ifSig],[0,1]);
            %end1 = strfind([ifSig, 0],[1,0]);
            sig(sig<=1e-20) = 1e-20;
            sig(sig>=0.05) = 1;
            logSig = log(sig)/log(1e-20);
            for ll=1:numel(sig)
                plot(t(ll)+dt*[-0.5 0.5], [104 104],'-','Color',[1 1-logSig(ll) 1-logSig(ll)],'LineWidth',5);
                %plot(t(ll)+dt*[-0.5 0.5], [95 95],'-','Color',[1 sig(ll) sig(ll)],'LineWidth',5);
            end

        else
            % control regression exist
            fieldName = ['trigEvent',num2str(currPredictor-1)];
            
            % find the continuous significant part
%             ifSig = (input_ctrl.(fieldName).boothigh < percentage);
%             start1 = strfind([0, ifSig'],[0,1]);
%             end1 = strfind([ifSig', 0],[1,0]);
% 
%             for ll = 1:length(start1)
%                 plot(t(start1(ll):end1(ll)), percentage(start1(ll):end1(ll)),'-','Color',[1 0.5 0],'LineWidth',5);
%             end
            sig = input_ctrl.(fieldName).pval;
            sig(sig<=1e-3) = 1e-3;
            %sig(sig>=0.05) = 1;
            logSig = log(sig)/log(1e-3);
            for ll=1:numel(sig)
                plot(t(ll)+dt*[-0.5 0.5], [95 95],'-','Color',[1 1-logSig(ll) 1-logSig(ll)],'LineWidth',5);
                %plot(t(ll)+dt*[-0.5 0.5], [95 95],'-','Color',[1 sig(ll) sig(ll)],'LineWidth',5);
            end

        end
    end
    
    if nInteraction > 0
        for l = 1:nInteraction
            for k=1:1+nback
                
                currPredictor=1+nPredictor*(1+nback)+(l-1)*(1+nback)+k;
                
                subplot(panelv,1+nback,currPredictor-1); hold on;
                %patch([t(1) t(end) t(end) t(1)],[0 0 100*pvalThresh 100*pvalThresh],[0.7 0.7 0.7],'EdgeColor','none');
                percentage = 100*sum(pval(:,currPredictor,:)<pvalThresh,3)/nCells;
                plot(t,percentage,'k.-','MarkerSize',15);
                plot([0 0],[0 100],'k','LineWidth',1);
                xlim([floor(t(1)) ceil(t(end))]);
                xticks([floor(t(1)):1:ceil(t(end))]);
                ylim([0 100]);
                title([tlabel{currPredictor-1}]);
                if if_xlabel == 1
                    xlabel(xtitle);
                end
                if if_ylabel == 1
                    ylabel('Fraction of sessions (%)');
                end

                
                %identifying significant points via binomial test
                if isempty(input_ctrl)
                    prop=[]; sig=[];
                    for ll=1:numel(t)
                        prop(ll)=sum(pval(ll,currPredictor,:)<pvalThresh,3)/nCells; %proportion of cells that is significant
                        [p]=myBinomTest(prop(ll)*nCells,nCells,pvalThresh);  %binomial test to see if this proportion is significantly different from the p-value threshold
                        sig(ll)=p;
                    end
                    for ll=1:numel(sig)
                        if sig(ll)<pvalThresh
                            plot(t(ll)+dt*[-0.5 0.5],[100 100],'k-','LineWidth',5);
                        end
                    end
                else
                    fieldName = ['trigEvent',num2str(currPredictor-1)];
                    for ll=1:numel(percentage)
                        if input_ctrl.(fieldName).boothigh(ll) < percentage(ll)
                            plot(t(ll)+dt*[-0.5 0.5],[100 100],'k-','LineWidth',5);
                        end
                    end
                end
            end
        end
        
    end
    
    %% for 4*3 (row*col)
    if nPredictor == 10 || nPredictor == 12
        pos = get(gca, 'Position');
        pos(4) = 0.13;
        pos(3) = 0.23;
        set(gca, 'Position', pos)
    elseif nPredictor == 14
    % for 4*4
        pos = get(gca, 'Position');
        pos(4) = 0.13;
        pos(3) = 0.15;
        set(gca, 'Position', pos)
    %% for 3*3
    elseif nPredictor == 9
        pos = get(gca, 'Position');
        pos(4) = 0.21;
        pos(3) = 0.23;
        set(gca, 'Position', pos)
    end
end
end