function xcorrData = getRegSelCorr(reg_cr,label,Ind, thresh, savefluofigpath)

% get the corresponding regression results
regCoeff = [];regPval = [];
reg_cr=reg_cr.reg_cr;
for ii = 1:length(reg_cr)
    regCoeff = [regCoeff;reg_cr{ii}.coeff(:,Ind)'];
    regPval = [regPval;reg_cr{ii}.pval(:,Ind)'];
end

Mat = selectRg2D(reg_cr,Ind,[]);

% take in the data, identify the significant time period

Mat.sig = Mat.pval;
Mat.sig(Mat.pval<thresh) = 1; Mat.sig(Mat.pval>=thresh) = 0;

% set not significant coefficient to NaN
Mat.sigCoeff = Mat.coeff;
Mat.sigCoeff(Mat.pval>=thresh) = 0;

% calculate frame-by-frame correlation, get the Index with maximum
% correlation achieved
sigCorr = zeros(size(Mat.sig,1)*2-1,size(Mat.sig,2)*2-1,size(Mat.sig,3)-1);
sigMaxInd = zeros(2,size(Mat.sig,3)-1);
sigMaxDis = zeros(1,size(Mat.sig,3)-1);  % find the coordinate of the maximum crosscorrelation, calculate its distance from the center
sigMaxValue = zeros(1,size(Mat.sig,3)-1);

centerCoor = [size(Mat.sig,1), size(Mat.sig,2)];

for ii = 1:size(Mat.sig,3)-1
    mat1 = Mat.sig(:,:,ii); mat2 = Mat.sig(:,:,ii+1);
    
    if all(mat1(:)==0) | all(mat2(:)==0)
        sigCorr(:,:,ii) = NaN;
        sigMaxInd(:,ii)=[NaN,NaN];
        sigMaxDis(ii) = NaN;
        sigMaxValue(ii) = NaN;
    else
        temp = normxcorr2(Mat.sig(:,:,ii),Mat.sig(:,:,ii+1));
        sigCorr(:,:,ii) = temp;
        [ssr,snd] = max(temp(:));
        sigMaxValue(ii) = ssr;
        [ij,ji] = ind2sub(size(temp),snd);
        sigMaxInd(:,ii)=[ij,ji];
        sigMaxDis(ii) = sqrt((ij-centerCoor(1))^2+(ji-centerCoor(2))^2);
    end
end

%% use significant coefficient only
coeffCorr = zeros(size(Mat.coeff,1)*2-1,size(Mat.coeff,2)*2-1,size(Mat.coeff,3)-1);
coeffMaxInd = zeros(2,size(Mat.coeff,3)-1);
coeffMaxDis = zeros(1,size(Mat.coeff,3)-1);  % find the coordinate of the maximum crosscorrelation, calculate its distance from the center
coeffMaxValue = zeros(1, size(Mat.coeff,3)-1);

for ii = 1:size(Mat.coeff,3)-1
    mat1 = Mat.sigCoeff(:,:,ii); mat2 = Mat.sigCoeff(:,:,ii+1);
    
    if all(mat1(:)==0) | all(mat2(:)==0)
        coeffCorr(:,:,ii) = NaN;
         coeffMaxInd(:,ii)=[NaN,NaN];
        coeffMaxDis(ii) = NaN;
        coeffMaxValue(ii) = NaN;
    else
        
        temp = normxcorr2(Mat.sigCoeff(:,:,ii),Mat.sigCoeff(:,:,ii+1));
        coeffCorr(:,:,ii) = temp;
        [ssr,snd] = max(temp(:));
        coeffMaxValue(ii) = ssr;
        [ij,ji] = ind2sub(size(temp),snd);
        coeffMaxInd(:,ii)=[ij,ji];
        coeffMaxDis(ii) = sqrt((ij-centerCoor(1))^2+(ji-centerCoor(2))^2);
    end
end

t= -2.95:0.1:4.85;

% plot the cross-correlation value of significance and coefficient
figure;
subplot(3,1,1)
plot(reg_cr{1}.regr_time,sum(regPval<thresh,1)/size(regPval,1),'-k');
ylabel('Significant Grids');
title(['Fraction of significant grids (',label,')']);
set(gca,'box','off')

subplot(3,1,2)
plot(t,squeeze(sigCorr(28,28,:)),'-r');
hold on;plot(t,squeeze(coeffCorr(28,28,:)),'-k');
set(gca,'box','off');
ylabel('Coefficient');
title(['Cross-correlation coefficient (',label,')']);
legend({'Significance','Coefficient'});
legend('box','off');
ylim([-0.2,1]);

subplot(3,1,3)
plot(t,sigMaxDis,'-r');
hold on;plot(t,coeffMaxDis,'-k');
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('Offset');
title(['Offset of maximum correlation coefficient (',label,')']);
legend({'Significance','Coefficient'});
legend('box','off');
print(gcf,'-dpng',fullfile(savefluofigpath,[label,' xcorrCoeff']));
saveas(gcf, fullfile(savefluofigpath,[label,' xcorrCoeff']), 'fig');

% plot the maximum cross correlation coefficient;
figure;
plot(t,sigMaxValue,'-r');
hold on;plot(t,coeffMaxValue,'-k');
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('Frame-by-frame cross-correlation');
title(['Max cross-correlation coefficient (',label,')']);
legend({'Significance','Coefficient'});
legend('box','off');
ylim([-0.2,1]);
print(gcf,'-dpng',fullfile(savefluofigpath,[label,'max cross-correaltion']));
saveas(gcf, fullfile(savefluofigpath,[label,'max cross-correaltion']), 'fig');

% plot the distance from center point of maximum crosscorrelation


% save the data into a structure
xcorrData.sigCorr = sigCorr;
xcorrData.sigMaxDis = sigMaxDis;
xcorrData.sigMaxInd = sigMaxInd;
xcorrData.sigMaxValue = sigMaxValue;
xcorrData.regPval = regPval;
xcorrData.regCoeff = regCoeff;
xcorrData.coeffCorr = coeffCorr;
xcorrData.coeffMaxDis = coeffMaxDis;
xcorrData.coeffMaxInd = coeffMaxInd;
xcorrData.coeffMaxValue = coeffMaxValue;

close all;
