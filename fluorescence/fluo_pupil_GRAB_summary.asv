function fluo_pupil_GRAB_summary(dataIndex)

% summarize the results of pupil/fluorescent correlation
% GRAB_Ach, before MP, 5,8,10,12

nFiles = size(dataIndex,1);
corrCoeffMat = [];
corrPMat = [];
for ii = 1:nFiles
     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    savematname = fullfile(fn_beh.folder,'analysis-pupil','corrPupFluo.mat');
    load(savematname);
    corrCoeffMat = [corrCoeffMat,corrCoeff(:)'];
    corrPMat = [corrPMat, corrP(:)'];
end

% plot the correlation histogram with significance indicated
maxCorr = max(corrCoeffMat);
minCorr = min(corrCoeffMat);
binSize = 0.01;
minBin = floor(minCorr/binSize)*binSize; maxBin = ceil(maxCorr/binSize)*binSize;
binCenter = minBin+binSize/2:binSize:maxBin-binSize/2;
numBin = length(binCenter);
binMat = zeros(numBin,1); binSig = zeros(numBin,1);

for cc = 1:length(corrCoeffMat)
    ind = floor((corrCoeffMat(cc)-minBin)/binSize)+1;
    binMat(ind) = binMat(ind)+1;
    if corrPMat(cc) < 0.05
        binSig(ind) = binSig(ind)+1;
    end
end

figure;bar(binCenter,binMat,'EdgeColor',[0.7,0.7,0.7],'FaceColor',[0.7,0.7,0.7])
hold on;bar(binCenter,binSig,'EdgeColor',[0.929,0.694,0.125],'FaceColor',[0.929,0.694,0.125]);
set(gca,'box','off');

