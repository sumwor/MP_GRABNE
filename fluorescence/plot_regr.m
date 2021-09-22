function plot_regr(input,pvalThresh,tlabel,xtitle)
% % plot_regr %
%PURPOSE:   Plot results from multiple linear regression
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   input:        Structure generated by linear_regr().
%   pvalThresh:   Threshold value to deem whether p-values are significant
%   tlabel:       Text to put as title of the plot.
%   xtitle:       Text to put as the label for x-axis.

%% setup
t=input{1}.regr_time;
dt=nanmean(diff(t));

nCells=numel(input);
for j=1:nCells
    pval(:,:,j)=input{j}.pval;
end

nPredictor=input{1}.numPredictor;

%number of row of panels
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

nback=input{1}.nback;

%% plot results

figure;

for l=1:nPredictor
    for k=1:1+nback
        currPredictor=1+(l-1)*(1+nback)+k; %first +1 because first term is bias
        
        subplot(panelv,1+nback,currPredictor-1); hold on;
        patch([t(1) t(end) t(end) t(1)],[0 0 100*pvalThresh 100*pvalThresh],[0.7 0.7 0.7],'EdgeColor','none');
        plot(t,100*sum(pval(:,currPredictor,:)<pvalThresh,3)/nCells,'k.-','MarkerSize',30);
        plot([0 0],[0 100],'k','LineWidth',1);
        xlim([t(1) t(end)]);
        ylim([0 40]);
        title(tlabel{currPredictor-1});
        if l == nPredictor && nInteraction == 0
            xlabel(xtitle);
        end
        if k == 1
            ylabel('Fraction of cells (%)');
        end
        
        %identifying significant points via binomial test
        sig=[];
        for ll=1:numel(t)
            [p]=myBinomTest(sum(pval(ll,currPredictor,:)<pvalThresh,3),nCells,pvalThresh);
            sig(ll)=p;
        end
        for ll=1:numel(sig)
            if sig(ll)<pvalThresh
                plot(t(ll)+dt*[-0.5 0.5],[35 35],'k-','LineWidth',5);
            end
        end
    end
end

% plot interaction terms, if present
if nInteraction > 0
    for l = 1:nInteraction
        for k=1:1+nback
            
            currPredictor=1+nPredictor*(1+nback)+(l-1)*(1+nback)+k;
            
            subplot(panelv,1+nback,currPredictor-1); hold on;
            patch([t(1) t(end) t(end) t(1)],[0 0 100*pvalThresh 100*pvalThresh],[0.7 0.7 0.7],'EdgeColor','none');
            plot(t,100*sum(pval(:,currPredictor,:)<pvalThresh,3)/nCells,'k.-','MarkerSize',30);
            plot([0 0],[0 100],'k','LineWidth',1);
            xlim([t(1) t(end)]);
            ylim([0 40]);
            title([tlabel{currPredictor-1}]);
            xlabel(xtitle);
            if k == 1
                ylabel('Fraction of cells (%)');
            end
            
            %identifying significant points via binomial test
            prop=[]; sig=[];
            for ll=1:numel(t)
                prop(ll)=sum(pval(ll,currPredictor,:)<pvalThresh,3)/nCells; %proportion of cells that is significant
                [p]=myBinomTest(prop(ll)*nCells,nCells,pvalThresh);  %binomial test to see if this proportion is significantly different from the p-value threshold
                sig(ll)=p;
            end
            for ll=1:numel(sig)
                if sig(ll)<pvalThresh
                    plot(t(ll)+dt*[-0.5 0.5],[35 35],'k-','LineWidth',5);
                end
            end
        end
    end
end

end