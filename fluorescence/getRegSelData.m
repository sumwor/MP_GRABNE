function RegData = getRegSelData(reg_cr,label,Ind, coeffThresh, sigThresh, savefluofigpath)

Mat = selectRg2D(reg_cr,Ind,[]);

% take in the data, identify the significant time period
nCells = size(Mat.pval,1)*size(Mat.pval,2);
fracSig = squeeze(sum(sum(Mat.pval<sigThresh.alpha,1),2)/nCells);
sig = zeros(1,numel(fracSig));

[p]=myBinomTest(squeeze(sum(sum(Mat.pval<sigThresh.alpha,1),2)),nCells,sigThresh.alpha);
ifSig = p<sigThresh.alpha;
%TF = ischange(double(ifSig));
ipt = findchangepts(double(ifSig),'MaxNumChanges',2);

% if only one change point found, determine the significant part
if length(ipt) == 1
    mean1 = mean(ifSig(1:ipt));mean2 = mean(ifSig(ipt:end));
    if mean2>mean1
        sigRange = [ipt:length(ifSig)];
    else
        sigRange = [1: ipt];
    end
elseif isempty(ipt)
    sigRange = [1:length(ifSig)];
else
    if mean(ifSig(ipt(1):ipt(2))) > mean(ifSig([1:ipt(1),ipt(2):end]))
        sigRange = [ipt(1):ipt(2)];  % need to check if it is true
    else
        sigRange = [1:ipt(1),ipt(2):length(ifSig)];
    end
end

%coefficient mask with significant grids only
RegData.coeffMask = getSignedSelMask(Mat,sigThresh.alpha,label, savefluofigpath);
[RegData.negMask,RegData.posMask] = getPNMask(RegData.coeffMask,coeffThresh);
isSig = 0;
% plotMaskBoundary(RegData.coeffMask,RegData.negMask,RegData.posMask,label,savefluofigpath,isSig);


% plot absolut coefficient with significant grids onlu
isSig = 2;
% get significant coefficient only
sigCoeff = Mat.coeff;
% determine the significant grids (threshold: at least 3 consecutive significant point
% or 10% of the time points are significant
sigMat = zeros(size(Mat.coeff,1),size(Mat.coeff,2));
for ss = 1:size(Mat.coeff,1)
    for uu = 1:size(Mat.coeff,2)
        sigVec = squeeze(Mat.pval(ss,uu,:)<sigThresh.alpha);
        start1 = strfind([0,sigVec'==1],[0 1]);
        end1 = strfind([sigVec'==1,0],[1 0]);
        if sum((end1-start1+1)>=3) | sum(sigVec) >= 8
            sigMat(ss,uu) = 1;
        end
    end
end

for zz = 1:size(Mat.coeff,3)
    temp = squeeze(Mat.coeff(:,:,zz));
    temp(~sigMat) = NaN;
    sigCoeff(:,:,zz) = temp;
end
% coeffBound.min = floor(min(min(min(sigCoeff)))*10)/10;
% coeffBound.max = ceil(max(max(max(sigCoeff)))*10)/10;
savepath = fullfile(savefluofigpath,'sigCoeff');
if ~exist(savepath)
    mkdir(savepath)
end
% plotMaskBoundary(sigCoeff,[],[],label,savepath,isSig);

% get the number of time points with pos/neg correlation coefficient

RegData.coeffSignSig = NaN(size(Mat.coeff,1),size(Mat.coeff,2),2);
t=-2.95:0.1:4.95;
for ss = 1:size(Mat.coeff,1)
    for uu = 1:size(Mat.coeff,2)
        if sigMat(ss,uu) == 1
%             figure;
%             plot(t,squeeze(sigCoeff(ss,uu,:)));
            
% only count significant time point?
            RegData.coeffSignSig(ss,uu,1) = sum(squeeze(sigCoeff(ss,uu,:)>0 & Mat.pval(ss,uu,:)<sigThresh.alpha))/sum(squeeze(Mat.pval(ss,uu,:)<sigThresh.alpha));
            RegData.coeffSignSig(ss,uu,2) = -sum(squeeze(sigCoeff(ss,uu,:)<0& Mat.pval(ss,uu,:)<sigThresh.alpha))/sum(squeeze(Mat.pval(ss,uu,:)<sigThresh.alpha));
        end
    end
end
RegData.signRel = RegData.coeffSignSig(:,:,1)./(-RegData.coeffSignSig(:,:,2));
% figure;histogram(RegData.signRel)

% significance mask
tlabel1 = [label,' significance'];
xtitle = 'Time from cue{s}';
colorRange=[0 0.01];
% t=-2.95:0.1:4.95;
% plot_signficance(Mat.pval,t, tlabel1,xtitle,colorRange,savefluofigpath)
%         
%         
% % plot average significance 
RegData.sigMask = getSigSelMask(Mat.pval, sigThresh.alpha, label,savefluofigpath);

[RegData.notSigMask] = getSigMask(RegData.sigMask,sigThresh.value);
% isSig=1;
% plotMaskBoundary(RegData.sigMask,RegData.notSigMask,[],label,savefluofigpath,isSig);

% make a video
colorRange = [-0.05,0.05];
video_coeff(sigCoeff,label,colorRange,savepath);