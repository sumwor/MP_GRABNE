function plot_gl_regr_fit(input_fit,input,params)
% % plot_gl_regr_fit %
%PURPOSE:   Plot predicted versus measured signal for GLM analysis
%AUTHORS:   AC Kwan 170901
%
%INPUT ARGUMENTS
%   input_fit:    Structure generated by gl_regr().
%   input:        Structure generated by gl_regr().
%   params:       Parameters for analysis and plotting.

figure;
subplot(2,1,1); hold on;
idx = ~isnan(input_fit.ytest);
plot(1:sum(idx),input_fit.ytest(idx),'k-','LineWidth',2);
plot(1:sum(idx),input_fit.ytest_fit(idx),'b-','LineWidth',2);
legend('Measured','Predicted','location','northeast');
ylabel('dF/F');
xlabel('Frame in test data set');
title(['Corr. coefficient = ' num2str(input.ytest_CC)]);

end