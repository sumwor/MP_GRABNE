function stats=algo_logreg_CRInt(stats,params,x)
% % algo_logreg_CRInt %
%PURPOSE:   Simulate a player based on logistic regression
%           regressors are choice and choice x reward
%AUTHORS:   AC Kwan 180517
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       bias - bias term
%       Ch - choice
%       RC - choice x reward
%   x:      which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

% make sure the params are column vectors
b_choice=params.Ch(:);
b_int=params.RC(:);

% how many trials to look back for logistic regression
step_back = numel(params.RC);

%probability of choosing left based on logistic regression
if stats.currTrial == 1   %if this is the first trial
    stats.pl(1,x) = 0.5;  %choose randomly
else
    if stats.currTrial>step_back
        backIdx = step_back;  %if at least number of trials permit, use step_back
    else
        backIdx = stats.currTrial-1;        %else use what is available
    end
    
    %choice and reward from the past trials under consideration
    c = stats.c(stats.currTrial-1:-1:stats.currTrial-backIdx,x);
    r = stats.r(stats.currTrial-1:-1:stats.currTrial-backIdx,x);
    
    % create regressor vectors
    Ch=zeros(length(c),1);   % right = 1; left = -1
    Ch=1*(c==1) + (-1)*(c==-1);

    Re=zeros(length(c),1);   % reward = 1; unrewarded = -1
    Re=r.*(r>0) + (-1)*(r==0);  %allows for reward size r!=1
%    Re=1*(r==1) + (-1)*(r==0); 
    
    RC=zeros(length(c),1);   % choice x reward interaction (right/reward=1; right/unreward=-1; left/reward=-1; left/unreward=1)
    RC=c.*r;
    
    b_temp = params.bias;    %bias term
    b_temp = b_temp + nansum(b_choice(1:backIdx).*Ch);    %choice terms
    b_temp = b_temp + nansum(b_int(1:backIdx).*RC);       %choice x reward terms
    stats.pl(stats.currTrial) = 1 - exp(b_temp)/(1+exp(b_temp));  %p(left) = 1 - p(right)
end

end
