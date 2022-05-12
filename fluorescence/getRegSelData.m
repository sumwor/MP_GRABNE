function RegData = getRegSelData(reg_cr,label,Ind, coeffThresh, sigThresh, savefluofigpath)

Mat = selectRg2D(reg_cr,Ind,[]);

% take in the data, identify the significant time period
% nCells = size(Mat.pval,1)*size(Mat.pval,2);
% fracSig = squeeze(sum(sum(Mat.pval<sigThresh.alpha,1),2)/nCells);
% sig = zeros(1,numel(fracSig));
% %for ll=1:numel(fracSig)
% [p]=myBinomTest(squeeze(sum(sum(Mat.pval<sigThresh.alpha,1),2)),nCells,sigThresh.alpha);
% ifSig = p<sigThresh.alpha;
% %TF = ischange(double(ifSig));
% ipt = findchangepts(double(ifSig),'MaxNumChanges',2);
% 
% % if only one change point found, determine the significant part
% if length(ipt) == 1
%     mean1 = mean(ifSig(1:ipt));mean2 = mean(ifSig(ipt:end));
%     if mean2>mean1
%         sigRange = [ipt:length(ifSig)];
%     else
%         sigRange = [1: ipt];
%     end
% elseif isempty(ipt)
%     sigRange = [1:length(ifSig)];
% else
%     if mean(ifSig(ipt(1):ipt(2))) > mean(ifSig([1:ipt(1),ipt(2):end]))
%         sigRange = [ipt(1):ipt(2)];  % need to check if it is true
%     else
%         sigRange = [1:ipt(1),ipt(2):length(ifSig)];
%     end
% end

% coefficient mask with significant grids only
RegData.coeffMask = getSignedSelMask(Mat,sigThresh.alpha,label, savefluofigpath);
[RegData.negMask,RegData.posMask] = getPNMask(RegData.coeffMask,coeffThresh);
isSig = 0;
plotMaskBoundary(RegData.coeffMask,RegData.negMask,RegData.posMask,label,savefluofigpath,isSig);


% plot absolut coefficient with significant grids onlu
isSig = 2;
% get significant coefficient only
sigCoeff = Mat.coeff;
sigCoeff(Mat.pval>=sigThresh.alpha) = 0;
coeffBound.min = floor(min(min(min(sigCoeff)))*10)/10;
coeffBound.max = ceil(max(max(max(sigCoeff)))*10)/10;
plotMaskBoundary(sigCoeff,[],[],label,savefluofigpath,isSig);


% significance mask
tlabel1 = [label,' significance'];
xtitle = 'Time from cue{s}';
colorRange=[0 0.01];
t=-2.95:0.1:4.95;
plot_signficance(Mat.pval,t, tlabel1,xtitle,colorRange,savefluofigpath)
        
        
% plot average significance 
RegData.sigMask = getSigSelMask(Mat.pval, sigThresh.alpha, label,savefluofigpath);

[RegData.notSigMask] = getSigMask(RegData.sigMask,sigThresh.value);
isSig=1;
plotMaskBoundary(RegData.sigMask,RegData.notSigMask,[],label,savefluofigpath,isSig);
