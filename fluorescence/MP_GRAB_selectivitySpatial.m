function MP_GRAB_selectivitySpatial(dataIndex);

% spatial properties of choice/outcome selectivity
%% probabily need more grids - pixel selectivity
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);
colorRange = [-1 1];

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
        
        load(fullfile(fn_fluo.folder,fn_fluo.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        saveSelName = fullfile(savematpath,'select_norm.mat');  % random forest
        %savePrevSelName = fullfile(savematpath,'prev_select_norm.mat'); 
        if exist(saveSelName) 
            load(saveSelName);
            
            % put the data in spatial organization
            edgelength = sqrt(numel(choicesel));
            choiceselMat = zeros(edgelength,edgelength,length(choicesel{1}.signal));
            
            
            for cc = 1:numel(choicesel)
                if mod(cc,edgelength) == 0
                    Ind2 = edgelength;
                else
                    Ind2 = mod(cc,(edgelength));
                end
                if mod(cc,edgelength) == 0
                    Ind1 = cc/edgelength;
                else
                    Ind1 = floor(cc/edgelength)+1;
                end
                
                % change left-right to contra-ipsi
              
                    if strcmp(dataIndex.RecordingSite{ii},'left')
                         choiceselMat(Ind1, Ind2,:) = -choicesel{cc}.signal;
                    else
                        choiceselMat(Ind1, Ind2,:) = choicesel{cc}.signal;
                    end
            end
            
            
            %% calculate the area with bwareaopen()
            posMask = zeros(edgelength,edgelength,size(choiceselMat,3));
            negMask = zeros(edgelength,edgelength,size(choiceselMat,3));
            mergeMask = zeros(edgelength,edgelength,size(choiceselMat,3));
            for tt = 1:size(choiceselMat,3)
                posBW = (choiceselMat(:,:,tt)>0); %Get logical mask of pixels exceeding threshold
                posBW = ~bwareaopen(~posBW,50,8); %Remove small holes from pixel mask
                posMask(:,:,tt) = bwareafilt(posBW,6,4);
                negBW = (choiceselMat(:,:,tt)<0); %Get logical mask of pixels exceeding threshold
                negBW = ~bwareaopen(~negBW,50,8); %Remove small holes from pixel mask
                negMask(:,:,tt) = -bwareafilt(negBW,6,4);

                % merge the pos/neg Mask (pos px = 1; neg px = -1);
                tempmerge = zeros(edgelength,edgelength);
                temppos = posMask(:,:,tt); tempneg = negMask(:,:,tt);
                tempmerge(temppos==1) =1;
                tempmerge(tempneg==-1) = -1;
                tempmerge(temppos==1 & tempneg==-1) = 0;
                mergeMask(:,:,tt) = tempmerge;   
            end

            figure;
            subplot(1,2,1)
            image(mean(posMask,3),'CDataMapping','scaled');
            axis square;
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]); 
            title('Average positive choice selecvitity mask');
            subplot(3,20,60);
            image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]);
            print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-pos mask'));
            saveas(gcf, fullfile(savefluofigpath,'c selecvitity-pos mask'), 'fig');

            figure;
            subplot(1,2,1)
            image(mean(negMask,3),'CDataMapping','scaled');
            axis square;
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]);
             title('Average negative choice selecvitity mask');
            subplot(3,20,60);
            image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]);
            print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-neg mask'));
            saveas(gcf, fullfile(savefluofigpath,'c selecvitity-neg mask'), 'fig');
            
             figure;
            subplot(1,2,1)
            image(mean(mergeMask,3),'CDataMapping','scaled');
            axis square;
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]); 
            title('Average merged choice selecvitity mask');
            subplot(3,20,60);
            image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]);
            print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-merge mask'));
            saveas(gcf, fullfile(savefluofigpath,'c selecvitity-merge mask'), 'fig');


%% directly averaging the sign of every grid
% pos = 1; neg = -1
% signMask = zeros(edgelength,edgelength,size(choiceselMat,3));
% signMask(choiceselMat>0) = 1;
% signMask(choiceselMat<0) = -1;
% signMask(choiceselMat == 0) = 0;
% 
% figure;
% subplot(1,2,1)
% image(mean(signMask,3),'CDataMapping','scaled');
% axis square;
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);
% title('Average positive choice selecvitity mask');
% subplot(3,20,60);
% image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);
% print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-sign mask'));
% saveas(gcf, fullfile(savefluofigpath,'c selecvitity-sign mask'), 'fig');

%% determine the area with watershed

% meanSignMask = mean(signMask,3);
% %tt = bwareaopen(meanSignMask,10,8);
% figure;imagesc(meanSignMask)
% absMask = abs(meanSignMask);
% figure;imagesc(absMask,'CDataMapping','scaled');
% axis square;
% colormap(colors);caxis([colorRange(1) colorRange(2)]);
% 
% 
% %figure;histogram(meanSignMask(:));
% absMasktt = imhmin(absMask,0.2);
% figure;imagesc(absMasktt,'CDataMapping','scaled');
% axis square;
% colormap(colors);caxis([colorRange(1) colorRange(2)]);
% 
% % convert the figure to binary
% level = graythresh(absMasktt);
% BW = im2bw(absMasktt,level);
% 
% BW = ~bwareaopen(~BW,10,8);
% BW = bwareaopen(BW,10,8);
% D = bwdist(~BW);
% %D(~BW) = -Inf;
% L = watershed(D);
% figure;image(label2rgb(L,'jet','w'))
% imshow(D,[])
% title('Distance Transform of Binary Image')
% I = imhmin(meanSignMask,0.9); %20 is the height threshold for suppressing shallow minima
% figure;image(I,'CDataMapping','scaled')
% axis square;
% colormap(colors);caxis([colorRange(1) colorRange(2)]);
% 
% L = watershed(D);
% figure;imagesc(L)

            %% get the areas, plot x/y arerage as a function of time
            % how to split the area? using bwareaopen&bwareafilt? manually?
            % threshold : 0-0.3; 0.3-0.7; 0.7-1.0
                area =  mean(mergeMask,3);
                BW1 = (area<-0.2); %Get logical mask of pixels exceeding threshold
                BW1 = ~bwareaopen(~BW1,10,4); %Remove small holes from pixel mask
                AMask1 = bwareafilt(BW1,6,4);
                L1 = bwlabel(AMask1);
                %function to get separate masks of objects in L
                Mask1 = getSepMask(L1,10);
                
                BW2 = (area>0.2); %Get logical mask of pixels exceeding threshold
                BW2 = ~bwareaopen(~BW2,10,4); %Remove small holes from pixel mask
                AMask2 = bwareafilt(BW2,6,4);
                %CC = bwconncomp(AMask2);
                L2 = bwlabel(AMask2);
                Mask2 = getSepMask(L2,10);
                
%                 BW3 = (area<0.65 & area>0); %Get logical mask of pixels exceeding threshold
%                 BW3 = ~bwareaopen(~BW3,10,8); %Remove small holes from pixel mask
%                 AMask3 = bwareafilt(BW3,6,4);
%                 %CC = bwconncomp(AMask2);
%                 L3 = bwlabel(AMask3);
%                 Mask3 = getSepMask(L3,10);
% %                 
%                 BW4 = (area>0.65); %Get logical mask of pixels exceeding threshold
%                 BW4 = ~bwareaopen(~BW4,10,8); %Remove small holes from pixel mask
%                 AMask4 = bwareafilt(BW4,6,4);
%                 %CC = bwconncomp(AMask2);
%                 L4 = bwlabel(AMask4);
%                 Mask4 = getSepMask(L4,10);
                
                
                %% go through the different areas, plot the average choice selevitity over time
                figure;
                subplot(1,2,1)
                image(mean(mergeMask,3),'CDataMapping','scaled');
                axis square;
                colormap(colors);
                caxis([colorRange(1) colorRange(2)]);
                title('Average choice selecvitity mask');
                hold on;
                %get a function to plot different Masks
                plotbound(Mask1);plotbound(Mask2);%plotbound(Mask3);plotbound(Mask4);
                subplot(3,20,60);
                image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
                colormap(colors);
                caxis([colorRange(1) colorRange(2)]);
                print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-sign mask-division'));
                saveas(gcf, fullfile(savefluofigpath,'c selecvitity-sign mask-division'), 'fig');
                
                %% save the mask information
                savemaskpath = fullfile(savematpath,'choiceselMask.mat');
                save(savemaskpath,'Mask1','Mask2')
%                 %
%                 aveChoiceSel1 = getAveSel(choiceselMat,Mask1);
%                 aveChoiceSel2 = getAveSel(choiceselMat,Mask2);
%                 aveChoiceSel3 = getAveSel(choiceselMat,Mask3);
%                 aveChoiceSel4 = getAveSel(choiceselMat,Mask4);
%                 t = -2.95:0.1:4.95;
%                 figure;
%                 plotAveSel(t,aveChoiceSel1,[0 0.4470 0.7410]);
%                 plotAveSel(t,aveChoiceSel2,[0.4660 0.6740 0.1880]);
%                 plotAveSel(t,aveChoiceSel3,[0.9290 0.6940 0.1250]);
%                 plotAveSel(t,aveChoiceSel4,[0.6350 0.0780 0.1840]);
%                 title('Average choice selectivity according to areas');
%                 xlabel('Time from cue(s)');
%                 ylabel('Average choice selectivity');
%                  print(gcf,'-dpng',fullfile(savefluofigpath,'ave choice sel by area'));
%                 saveas(gcf, fullfile(savefluofigpath,'ave choice sel by area'), 'fig');
%                 
                    close all;
%% monte carlo simulation of
                % find the center of neg/pos areas
                % model the area by average density distribution, with
                % error
                % calculate the distribution of center distance
                % generate infinite plane by center distance and
                % distribution
               % poisson point process: number of points follows poisson
               % distribution
        end
    end
end
end