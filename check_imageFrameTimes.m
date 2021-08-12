function check_imageFrameTimes ( cells, stackInfo, savefluofigpath)
% % check_imageFrameTimes %
%PURPOSE:   Check the synchronization of imaging + behavioral data
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   cells:          Structure generated by calc_dFF().
%   stackInfo:      content of the stackInfo.mat file for this imaging session
%

%% Make sure the image frame times are reasonable
% Time between frames should be around dt

tlabel=char(stackInfo.savFile_name);
dt=1/stackInfo.frameRate;

figure;
subplot(2,1,1); hold on;
plot(1:numel(cells.t)-1, diff(cells.t),'k');
plot([1 numel(cells.t)-1], dt*[1.05 1.05],'r');
plot([1 numel(cells.t)-1], dt*[0.95 0.95],'r');
xlabel('Frames');
ylabel('Time between frames (s)');
title({tlabel;'Value should not deviate from +/-5% of nominal frame rate (red lines)'});

savefilename = fullfile(savefluofigpath, 'check-frametiming');
print(gcf,'-dpng',savefilename);    %png format
saveas(gcf, savefilename, 'fig');


end