function errorshade(x,h,l,c)
% % errorshade %
%PURPOSE:   Plot shaded area, e.g. for confidence intervals
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   x:  the independent variable
%   h:  the upper bound for the dependent variable
%   l:  the lower bound for the dependent variable
%   c:  the color of the shading, e.g., [0.7 0.7 0.7] for gray

%%
x=x(:)';
h=h(:)';
l=l(:)';

h=fill([x fliplr(x)],[h fliplr(l)],c,'LineStyle','none');