% GRABNE_patch

% find the spatial and temporal structure of adrenergic firing

filepath = "F:\GRABNE\tone\tone-pupil-1\mat";

cd(filepath);
matfiles = dir('*.mat');

IM = load(fullfile(matfiles(1).folder,matfiles(1).name));
figure;
for i = 1:912
    imagesc(IM.stack(:,:,i))
    axis equal
    pause(.002);
end

% pixel-wise cross correlation?
% binned cross0-correlation?

% use desampled frames first
desampledPath = 'F:\GRABNE\tone\tone-pupil-1';
desampledFile = 'file_0_DS8.tif';
desampledmat = loadtiffseq(desampledPath,desampledFile);

% frequency of the desampled video is ~ 2Hz
% binned cross-correlation?
meanInt = mean(desampledmat,3);
figure;imagesc(meanInt);axis equal;


binSize = 10; % 10 frames, ~ 5s long
nBin = floor(size(desampledmat,3)/binSize);
range = 100; % check -50 to + 50 pixels
% try pixal (380,335)
corrPix = zeros(range+1,range+1,nBin);
cX = 335; cY = 380;
for xx = cX-range/2:cX+range/2
    for yy = cY-range/2:cY+range/2
        for bb = 1:nBin
            corrMat = corrcoef(double(desampledmat(cX,cY,binSize*(bb-1)+1:binSize*bb)),double(desampledmat(xx,yy,binSize*(bb-1)+1:binSize*bb))); 
            corrPix(xx-cX+range/2+1,yy-cY+range/2+1,bb) = corrMat(1,2);
        end
    end
end

% exclude surronding pixels (16 pixel)
pxExclude = 16;
corrThresh = 0.2;
minPx = 300;  % ROIs less than 300 are not considered
% also don't go through every pixel. Instead do 1:10:end
iterStep = 10;
xIter = pxExclude+1:iterStep:size(desampledmat,1)-pxExclude; yIter = pxExclude+1:iterStep:size(desampledmat,1)-pxExclude;
cellMask = [];
tic
for xx = xIter
%     if mod(xx,10) == 0
%         display(xx)
%     end
    for yy = yIter
        
        % check if given pixel already in the BW 
        if (~isempty(cellMask) & sum(cellMask(xx,yy,:)) == 0) | isempty(cellMask)
            corrPixAll = zeros(size(desampledmat,1)-pxExclude*2,  size(desampledmat,2)-pxExclude*2);
            for xX = pxExclude+1:size(desampledmat,1)-pxExclude
                for yY = pxExclude+1:size(desampledmat,1)-pxExclude
                    corrMat = corrcoef(double(desampledmat(xX,yY,:)),double(desampledmat(xx,yy,:))); 
                    corrPixAll(xX,yY) = corrMat(1,2);
                end
            end
            BW = (corrPixAll>corrThresh); %Get logical mask of pixels exceeding threshold
            BW = ~bwareaopen(~BW,10,8);
            roiMask = bwareafilt(BW,1,4);
            if sum(sum(roiMask==1)) >= minPx
                cellMask = cat(3,cellMask, roiMask); %Keep only the largest region; 'conn'=4
            end
        end
    end
end
toc

corrPixAll = zeros(size(desampledmat,1)-pxExclude*2,  size(desampledmat,2)-pxExclude*2);
        for xx = 1:xIter
            for yy = 1:yIter
                corrMat = corrcoef(double(desampledmat(xx,yy,:)),double(desampledmat(cX,cY,:))); 
                corrPixAll(xx,yy) = corrMat(1,2);
            end
        end
        BW = (corrPixAll>corrThresh); %Get logical mask of pixels exceeding threshold
        BW = ~bwareaopen(~BW,10,8); 
        
        cellMask{xx,yy} = bwareafilt(BW,1,4);
        
 figure;imagesc(corrPixAll);axis equal
 colorbar
 figure;imagesc(BW);axis equal
 
 figure;histogram(corrPixAll)