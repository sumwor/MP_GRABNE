function [ output ] = linear_regr_lick ( signal, t, event, eventTime, trialSubset, params)
% % linear_regr %
%% change the regression equation to Sul et al.2011
%PURPOSE:   Multiple lienar regression
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   signal:         time-series signal (e.g., dF/F)
%   t:              time points corresponding to the signal
%   event:          event, dummy-coded (e.g., for choice, left=-1, right=1)
%                   %currently can handle up to 2 types of events, e.g., choice and outcome
%   eventTime:      the event times
%   trialSubset:    the subset of trials to investigate, all else set to NaN
%   params.window:  the time window around which to align signal
%   params.nback:   consider events from up to this # trials back


%OUTPUT ARGUMENTS
%   output:         structure containing regression analysis results
%
% To plot the output, use MP_plot_regr().

%% interpolate signal

window=params.window;

% interpolate the signal to a finer time scale
interdt=0.01;
intert=[t(1):interdt:t(end)]';

% t might contain similar time point
intersig=interp1(t,signal,intert);

% align signal to the event
% use window slightly wider than the regression, so regression analysis
% won't run into the boundaries of this variable
[sigbyTrial, tbyTrial]=align_signal(intert,intersig,eventTime,[window(1)-1 window(end)+1]);

%% the factors in the linear regression

output.numPredictor=size(event,2)+1;  %number of predictor
output.nback = params.nback;
output.interaction = params.interaction;

factors=[];
%for j=1:output.numPredictor-2
for j=1:output.numPredictor-1
    event(~trialSubset,j)=nan;        %if event is not part of trialSubset to analyze, set to NaN
    factors = [factors event(:,j)];
    
    if output.numPredictor == 1 || output.numPredictor == 2
        for k=1:output.nback
            event_kback=[NaN(k,1); event(1:end-k,j)];     %event k-back
            factors = [factors event_kback];
        end
    end
end




if output.numPredictor==1
    terms=[zeros(1,1+output.nback) ; eye(1+output.nback)]; %bias, c(n), c(n-1), c(n-2), c(n-3), ...
elseif output.numPredictor==2
    if params.interaction == true
        %bias, c(n), c(n-1), c(n-2), ..., r(n), r(n-1), r(n-2), ...,
        terms=[zeros(1,(1+output.nback)*output.numPredictor) ; eye((1+output.nback)*output.numPredictor)];
        %add interaction terms: c(n)xr(n), c(n-1)xr(n-1), ...
        terms=[terms; [eye(1+output.nback) eye(1+output.nback)]];
    else
        %bias, c(n), c(n-1), c(n-2), ..., r(n), r(n-1), r(n-2), ...
        terms=[zeros(1,(1+output.nback)*output.numPredictor) ; eye((1+output.nback)*output.numPredictor)];
    end
elseif output.numPredictor > 2 % this is the linear regression for RL latent variables
    %disp('Linear regression model');
    terms = [zeros(1, output.numPredictor);eye(output.numPredictor,output.numPredictor)];
end


%% parameters of the regression
step_dur = nanmean(diff(window));
step_size = step_dur;       %if step size is step duration, then doing this in non-overlapping windows
output.regr_time=[window(1)+step_dur/2:step_size:window(end)-step_dur/2]';

%% multiple linear regression with moving window
% to see which behavioral events explain signal variability

%% test
nLick = zeros(length(eventTime),numel(output.regr_time));
%nIpsiLick = zeros(length(eventTime),numel(output.regr_time));
%nContraLick = zeros(length(eventTime),numel(output.regr_time));
for jj=1:numel(output.regr_time)
    
    
    for tt = 1:length(eventTime)
        startTime = output.regr_time(jj)-step_dur/2 + params.trigTime(tt);
        endTime = output.regr_time(jj)+step_dur/2 + params.trigTime(tt);
        nLick(tt,jj) = sum(params.lick>=startTime & params.lick<endTime);
        
        %nIpsiLick(tt,jj) = sum(params.ipsilick>=startTime & params.ipsilick<endTime);
        %nContraLick(tt,jj) = sum(params.contralick>=startTime & params.contralick<endTime);
    end
end

warning('off', 'all');
warning('off','MATLAB:singularMatrix');
warning('off','stats:pvaluedw:ExactUnavailable');
for jj=1:numel(output.regr_time)
    idx1=sum(tbyTrial<=output.regr_time(jj)-step_dur/2);     %signal corresponding to current time step
    idx2=sum(tbyTrial<(output.regr_time(jj)+step_dur/2));
    tempsig=squeeze(nanmean(sigbyTrial(idx1:idx2,:),1));
    
    % add lick to the terms, count the time within time interval
        if sum(nLick(:,jj))>0
            temp_factors = [factors, nLick(:,jj)];
            temp_terms = terms;
            ifLick = 1;
        else
            temp_factors = factors;
            temp_terms = terms(1:end-1,1:end-1);
            ifLick = 0;
        end
%     if sum(nIpsiLick(:,jj))>0 && sum(nContraLick(:,jj))>0
%         temp_factors = [factors, nContraLick(:,jj), nIpsiLick(:,jj)];
%         temp_terms = terms;
%         lickEmpty = '' ; % make a record on which lick is absent
%     elseif sum(nIpsiLick(:,jj))>0 && sum(nContraLick(:,jj))==0
%         temp_factors = [factors, nIpsiLick(:,jj)];
%         temp_terms = terms(1:end-1,1:end-1);
%         lickEmpty = 'contra';
%     elseif sum(nIpsiLick(:,jj))==0 && sum(nContraLick(:,jj))>0
%         temp_factors = [factors, nContraLick(:,jj)];
%         temp_terms = terms(1:end-1,1:end-1);
%         lickEmpty = 'ipsi';
%     else
%         temp_factors = factors;
%         temp_terms = terms(1:end-2,1:end-2);
%         lickEmpty = 'both';
%     end
    
    
    %try
    warning('off', 'MATLAB:rankDeficientMatrix');
    stats=regstats(tempsig',temp_factors,temp_terms);
    mdl = fitlm(temp_factors,tempsig',temp_terms);
    tbl=anova(mdl);
    %for kk=1:output.numPredictor+1
    for kk=1:output.numPredictor+1
        
        if sum(nLick(:,jj))>0
            output.coeff(jj,kk)=stats.beta(kk);        %coefficient
            output.pval(jj,kk)=stats.tstat.pval(kk);
        else
            if kk < 16
                output.coeff(jj,kk)=stats.beta(kk);        %coefficient
                output.pval(jj,kk)=stats.tstat.pval(kk);
            else
                output.pval(jj,kk)=NaN;
                output.coeff(jj,kk)=NaN;
            end
        end
            
%         if isempty(lickEmpty)
%             output.coeff(jj,kk)=stats.beta(kk);        %coefficient
%             output.pval(jj,kk)=stats.tstat.pval(kk); %pvalue for coefficient
%         elseif strcmp(lickEmpty,'ipsi')
%             if kk < 17
%                 output.coeff(jj,kk)=stats.beta(kk);        %coefficient
%                 output.pval(jj,kk)=stats.tstat.pval(kk);
%             else
%                 output.coeff(jj,kk)=NaN;        %coefficient
%                 output.pval(jj,kk)=NaN;
%             end
%         elseif strcmp(lickEmpty,'contra')
%             if kk < 16
%                 output.coeff(jj,kk)=stats.beta(kk);        %coefficient
%                 output.pval(jj,kk)=stats.tstat.pval(kk);
%             elseif kk==16
%                 output.coeff(jj,kk)=NaN;        %coefficient
%                 output.pval(jj,kk)=NaN;
%             else
%                 output.coeff(jj,kk)=stats.beta(end);        %coefficient
%                 output.pval(jj,kk)=stats.tstat.pval(end);
%             end
%         elseif strcmp(lickEmpty,'both')
%             if kk < 16
%                 output.coeff(jj,kk)=stats.beta(kk);        %coefficient
%                 output.pval(jj,kk)=stats.tstat.pval(kk);
%             else
%                 output.coeff(jj,kk)=NaN;        %coefficient
%                 output.pval(jj,kk)=NaN;
%             end
%         end
        output.rsquare(jj,kk) = stats.rsquare;
        output.adjrsquare(jj,kk) = stats.adjrsquare;
        %             if kk == 1
        %                 output.SumSq(jj,kk) = tbl.SumSq(end);  % last one is error term
        %             else
        %                 output.SumSq(jj,kk) = tbl.SumSq(kk-1);
        %             end
    end
    output.mse(jj) = stats.mse;

    %     catch
    %         for kk=1:size(terms,1)
    %             output.coeff(jj,kk)=NaN;
    %             output.pval(jj,kk)=NaN;
    %             output.SumSq(jj,kk) = NaN;
    %             output.rsquare(jj,kk) = NaN;
    %             output.adjrsquare(jj,kk) = NaN;
    %         end
    %         output.mse(jj) = NaN;
    %     end
end
warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
