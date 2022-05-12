function plot_xcorrCoeff(input1,input2,input3,thresh,save_path_fluo,label)

input2 = getBootstrp(input2, 0, 0.05);
input3 = getBootstrp(input3, 0, 0.05);

figure;
subplot(3,1,1)
plot(input1.t,sum(sum(input1.coeff<thresh,1),3)/(size(input1.coeff,1)*size(input1.coeff,3)),'-k');
ylabel('Significant Grids');
title(['Fraction of significant grids (',label,')']);
set(gca,'box','off')

subplot(3,1,2)
plot(input2.t,input2.coeff_bootave);
hold on;
gray=[0.7 0.7 0.7];
errorshade(input2.t,input2.bootlow,input2.boothigh,gray);
hold on;
plot(input2.t,input2.coeff_bootave,'k-');
set(gca,'box','off');
ylabel('Coefficient');
title(['Cross-correlation coefficient (',label,')']);
%ylim([-0.2,1]);

subplot(3,1,3)
plot(input3.t,input3.coeff_bootave);
hold on;
gray=[0.7 0.7 0.7];
errorshade(input3.t,input3.bootlow,input3.boothigh,gray);
hold on;
plot(input3.t,input3.coeff_bootave,'k-');
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('Offset');
title(['Offset of maximum correlation coefficient (',label,')']);


print(gcf,'-dpng',fullfile(save_path_fluo,[label,' xcorrCoeff']));
saveas(gcf, fullfile(save_path_fluo,[label,' xcorrCoeff']), 'fig');
%bootstrap the results and plot it



% figure;
% for ii=1:size(input.coeff,1)
%     hold on;plot(input.t,input.coeff(ii,:));
% end