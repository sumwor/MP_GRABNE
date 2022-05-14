function MP_GRAB_selectivitySpatial(dataIndex)

% spatial properties of choice/outcome selectivity
%% probabily need more grids - pixel selectivity


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

%         saveSelName = fullfile(savematpath,'select_norm.mat');  % random forest
%         %savePrevSelName = fullfile(savematpath,'prev_select_norm.mat'); 
%         if exist(saveSelName) 
%             load(saveSelName);
%             
%             % put the data in spatial organization
%             edgelength = sqrt(numel(choicesel));
%             choiceselMat = zeros(edgelength,edgelength,length(choicesel{1}.signal));
%             
%             
%             for cc = 1:numel(choicesel)
%                 if mod(cc,edgelength) == 0
%                     Ind2 = edgelength;
%                 else
%                     Ind2 = mod(cc,(edgelength));
%                 end
%                 if mod(cc,edgelength) == 0
%                     Ind1 = cc/edgelength;
%                 else
%                     Ind1 = floor(cc/edgelength)+1;
%                 end
%                 
%                 % change left-right to contra-ipsi
%               
%                     if strcmp(dataIndex.RecordingSite{ii},'left')
%                          choiceselMat(Ind1, Ind2,:) = -choicesel{cc}.signal;
%                     else
%                         choiceselMat(Ind1, Ind2,:) = choicesel{cc}.signal;
%                     end
%             end
%             
%             
%             %% calculate the area with bwareaopen()
%             posMask = zeros(edgelength,edgelength,size(choiceselMat,3));
%             negMask = zeros(edgelength,edgelength,size(choiceselMat,3));
%             mergeMask = zeros(edgelength,edgelength,size(choiceselMat,3));
%             for tt = 1:size(choiceselMat,3)
%                 posBW = (choiceselMat(:,:,tt)>0); %Get logical mask of pixels exceeding threshold
%                 posBW = ~bwareaopen(~posBW,50,8); %Remove small holes from pixel mask
%                 posMask(:,:,tt) = bwareafilt(posBW,6,4);
%                 negBW = (choiceselMat(:,:,tt)<0); %Get logical mask of pixels exceeding threshold
%                 negBW = ~bwareaopen(~negBW,50,8); %Remove small holes from pixel mask
%                 negMask(:,:,tt) = -bwareafilt(negBW,6,4);
% 
%                 % merge the pos/neg Mask (pos px = 1; neg px = -1);
%                 tempmerge = zeros(edgelength,edgelength);
%                 temppos = posMask(:,:,tt); tempneg = negMask(:,:,tt);
%                 tempmerge(temppos==1) =1;
%                 tempmerge(tempneg==-1) = -1;
%                 tempmerge(temppos==1 & tempneg==-1) = 0;
%                 mergeMask(:,:,tt) = tempmerge;   
%             end
% 
%             figure;
%             subplot(1,2,1)
%             image(mean(posMask,3),'CDataMapping','scaled');
%             axis square;
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]); 
%             title('Average positive choice selecvitity mask');
%             subplot(3,20,60);
%             image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%             print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-pos mask'));
%             saveas(gcf, fullfile(savefluofigpath,'c selecvitity-pos mask'), 'fig');
% 
%             figure;
%             subplot(1,2,1)
%             image(mean(negMask,3),'CDataMapping','scaled');
%             axis square;
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%              title('Average negative choice selecvitity mask');
%             subplot(3,20,60);
%             image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%             print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-neg mask'));
%             saveas(gcf, fullfile(savefluofigpath,'c selecvitity-neg mask'), 'fig');
%             
%              figure;
%             subplot(1,2,1)
%             image(mean(mergeMask,3),'CDataMapping','scaled');
%             axis square;
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]); 
%             title('Average merged choice selecvitity mask');
%             subplot(3,20,60);
%             image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%             print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-merge mask'));
%             saveas(gcf, fullfile(savefluofigpath,'c selecvitity-merge mask'), 'fig');
% 
% 
% %% directly averaging the sign of every grid
% % pos = 1; neg = -1
% % signMask = zeros(edgelength,edgelength,size(choiceselMat,3));
% % signMask(choiceselMat>0) = 1;
% % signMask(choiceselMat<0) = -1;
% % signMask(choiceselMat == 0) = 0;
% % 
% % figure;
% % subplot(1,2,1)
% % image(mean(signMask,3),'CDataMapping','scaled');
% % axis square;
% % colormap(colors);
% % caxis([colorRange(1) colorRange(2)]);
% % title('Average positive choice selecvitity mask');
% % subplot(3,20,60);
% % image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
% % colormap(colors);
% % caxis([colorRange(1) colorRange(2)]);
% % print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-sign mask'));
% % saveas(gcf, fullfile(savefluofigpath,'c selecvitity-sign mask'), 'fig');
% 
% %% determine the area with watershed
% 
% % meanSignMask = mean(signMask,3);
% % %tt = bwareaopen(meanSignMask,10,8);
% % figure;imagesc(meanSignMask)
% % absMask = abs(meanSignMask);
% % figure;imagesc(absMask,'CDataMapping','scaled');
% % axis square;
% % colormap(colors);caxis([colorRange(1) colorRange(2)]);
% % 
% % 
% % %figure;histogram(meanSignMask(:));
% % absMasktt = imhmin(absMask,0.2);
% % figure;imagesc(absMasktt,'CDataMapping','scaled');
% % axis square;
% % colormap(colors);caxis([colorRange(1) colorRange(2)]);
% % 
% % % convert the figure to binary
% % level = graythresh(absMasktt);
% % BW = im2bw(absMasktt,level);
% % 
% % BW = ~bwareaopen(~BW,10,8);
% % BW = bwareaopen(BW,10,8);
% % D = bwdist(~BW);
% % %D(~BW) = -Inf;
% % L = watershed(D);
% % figure;image(label2rgb(L,'jet','w'))
% % imshow(D,[])
% % title('Distance Transform of Binary Image')
% % I = imhmin(meanSignMask,0.9); %20 is the height threshold for suppressing shallow minima
% % figure;image(I,'CDataMapping','scaled')
% % axis square;
% % colormap(colors);caxis([colorRange(1) colorRange(2)]);
% % 
% % L = watershed(D);
% % figure;imagesc(L)
% 
%             %% get the areas, plot x/y arerage as a function of time
%             % how to split the area? using bwareaopen&bwareafilt? manually?
%             % threshold : 0-0.3; 0.3-0.7; 0.7-1.0
%                 area =  mean(mergeMask,3);
%                 BW1 = (area<-0.2); %Get logical mask of pixels exceeding threshold
%                 BW1 = ~bwareaopen(~BW1,10,4); %Remove small holes from pixel mask
%                 AMask1 = bwareafilt(BW1,6,4);
%                 L1 = bwlabel(AMask1);
%                 %function to get separate masks of objects in L
%                 Mask1 = getSepMask(L1,10);
%                 
%                 BW2 = (area>0.2); %Get logical mask of pixels exceeding threshold
%                 BW2 = ~bwareaopen(~BW2,10,4); %Remove small holes from pixel mask
%                 AMask2 = bwareafilt(BW2,6,4);
%                 %CC = bwconncomp(AMask2);
%                 L2 = bwlabel(AMask2);
%                 Mask2 = getSepMask(L2,10);
%                 
% %                 BW3 = (area<0.65 & area>0); %Get logical mask of pixels exceeding threshold
% %                 BW3 = ~bwareaopen(~BW3,10,8); %Remove small holes from pixel mask
% %                 AMask3 = bwareafilt(BW3,6,4);
% %                 %CC = bwconncomp(AMask2);
% %                 L3 = bwlabel(AMask3);
% %                 Mask3 = getSepMask(L3,10);
% % %                 
% %                 BW4 = (area>0.65); %Get logical mask of pixels exceeding threshold
% %                 BW4 = ~bwareaopen(~BW4,10,8); %Remove small holes from pixel mask
% %                 AMask4 = bwareafilt(BW4,6,4);
% %                 %CC = bwconncomp(AMask2);
% %                 L4 = bwlabel(AMask4);
% %                 Mask4 = getSepMask(L4,10);
%                 
%                 
%                 %% go through the different areas, plot the average choice selevitity over time
%                 figure;
%                 subplot(1,2,1)
%                 image(mean(mergeMask,3),'CDataMapping','scaled');
%                 axis square;
%                 colormap(colors);
%                 caxis([colorRange(1) colorRange(2)]);
%                 title('Average choice selecvitity mask');
%                 hold on;
%                 %get a function to plot different Masks
%                 plotbound(Mask1);plotbound(Mask2);%plotbound(Mask3);plotbound(Mask4);
%                 subplot(3,20,60);
%                 image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%                 colormap(colors);
%                 caxis([colorRange(1) colorRange(2)]);
%                 print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-sign mask-division'));
%                 saveas(gcf, fullfile(savefluofigpath,'c selecvitity-sign mask-division'), 'fig');
%                 
%                 %% save the mask information
%                 savemaskpath = fullfile(savematpath,'choiceselMask.mat');
%                 save(savemaskpath,'Mask1','Mask2')
% %                 %
% %                 aveChoiceSel1 = getAveSel(choiceselMat,Mask1);
% %                 aveChoiceSel2 = getAveSel(choiceselMat,Mask2);
% %                 aveChoiceSel3 = getAveSel(choiceselMat,Mask3);
% %                 aveChoiceSel4 = getAveSel(choiceselMat,Mask4);
% %                 t = -2.95:0.1:4.95;
% %                 figure;
% %                 plotAveSel(t,aveChoiceSel1,[0 0.4470 0.7410]);
% %                 plotAveSel(t,aveChoiceSel2,[0.4660 0.6740 0.1880]);
% %                 plotAveSel(t,aveChoiceSel3,[0.9290 0.6940 0.1250]);
% %                 plotAveSel(t,aveChoiceSel4,[0.6350 0.0780 0.1840]);
% %                 title('Average choice selectivity according to areas');
% %                 xlabel('Time from cue(s)');
% %                 ylabel('Average choice selectivity');
% %                  print(gcf,'-dpng',fullfile(savefluofigpath,'ave choice sel by area'));
% %                 saveas(gcf, fullfile(savefluofigpath,'ave choice sel by area'), 'fig');
% %                 
%                     close all;
%% monte carlo simulation of
                % find the center of neg/pos areas
                % model the area by average density distribution, with
                % error
                % calculate the distribution of center distance
                % generate infinite plane by center distance and
                % distribution
               % poisson point process: number of points follows poisson
               % distribution

        
        warning('off','all')
        
           %% set threshold
        Thresh.pos = 0.2; Thresh.neg = -0.2; % threshold to find area with pos/neg modulation
        sigThresh.value = 0.05; sigThresh.alpha = 0.01;  % threshold to find area with significant modulation
        
        tic
        %% choice/outcome selectivity
%         saveSelName = fullfile(savematpath,'select_norm.mat');
%         %savePrevSelName = fullfile(savematpath,'prev_select_norm.mat');
%         savemaskpath = fullfile(savematpath,'selectivityMask.mat');  % file path to save the results
%         if exist(saveSelName)
%             load(saveSelName);
%             
%             if exist(savemaskpath) % if mask file exist
%                 load(savemaskpath)
%             end
%             % get data from selectivity
%             % put the data in spatial organization
%             if exist('choicesel','var') & ~exist('choiceselData', 'var') % if choice selectivity mask not computed
%                 choiceselMat = select2D(choicesel,dataIndex.RecordingSite{ii});
%                 % get signed average selectivity mask
%                 label = 'choice';
%                 choiceselData.mergeMask = getSignedSelMask(choiceselMat, label, savefluofigpath);
%                 % get the areas
%                 [choiceselData.negMask,choiceselData.posMask] = getPNMask(choiceselData.mergeMask, Thresh);
%                 plotMaskBoundary(choiceselData.mergeMask,choiceselData.negMask,choiceselData.posMask,label,savefluofigpath);
%                 % save the result
%             end
%             
%             if exist('outcomesel','var') & ~exist('outcomeselData', 'var') % if outcome selectivity mask not computed
%                 % outcome
%                 outcomeselMat = select2D(outcomesel,[]);
%                 
%                 % get signed average selectivity mask
%                 label = 'reward';
%                 outcomeselData.mergeMask = getSignedSelMask(outcomeselMat, label, savefluofigpath);
%                 % get the areas
%                 [outcomeselData.negMask,outcomeselData.posMask] = getPNMask(outcomeselData.mergeMask,Thresh);
%                 plotMaskBoundary(outcomeselData.mergeMask,outcomeselData.negMask,outcomeselData.posMask,label, savefluofigpath);
%                 
%             end
%             
%             save(savemaskpath,'choiceselData','outcomeselData')
% 
%         end
%         
        
        %% check other variables using regression results (should focus on significance rather than coefficient)
        
        % variable list:
        % regression 1: cn,rn,cn+1,cn-1,rn-1,cn*rn,average_r,cumulative_r
        % regression 2: dQ,dK,chosenQ,chosenK
        % regression 3: RPE, CKE,
        saveregmaskpath = fullfile(savematpath,'regressionMask.mat');  % file path to save the results
        
        if exist(saveregmaskpath)
            load(saveregmaskpath);
        end
        
        % need to get the boundaries of possible coefficient value of all
        % regressions
        saveCRName = fullfile(savematpath,'regCR_norm.mat');
         if exist(saveCRName)
            reg1 = load(saveCRName);
         end
        saveRLName = fullfile(savematpath,'regRL_norm.mat');
          if exist(saveRLName)
            reg2 = load(saveRLName);
         end
        saveRPEName = fullfile(savematpath,'regRL_RPE_norm.mat');
          if exist(saveRPEName)
            reg3 = load(saveRPEName);
          end
         
          % get the boundary of coefficient
%           maxCoeff = 0;minCoeff = 0;
%           for rr = 1:length(reg1.reg_cr)
%               tempMax = max([max(reg1.reg_cr{rr}.coeff(:)),max(reg2.reg_cr{rr}.coeff(:)),max(reg3.reg_cr{rr}.coeff(:))]);
%               tempMin = min([min(reg1.reg_cr{rr}.coeff(:)),min(reg2.reg_cr{rr}.coeff(:)),min(reg3.reg_cr{rr}.coeff(:))]);
%               if tempMax > maxCoeff
%                   maxCoeff = tempMax;
%               end
%               if tempMin < minCoeff
%                   minCoeff = tempMin;
%               end
%           end
%                       coeffBound.min = floor(minCoeff*10)/10;
%                       coeffBound.max = ceil(maxCoeff*20)/10;
                          
        %% regression 1 -------------------------------------------------------------------------
             % get the choice from linear regression
        if ~exist('choiceRegData','var') % if choice regression mask not computed
            label = 'choice';
            choiceInd = 3;

            choiceRegData = getRegSelData(reg1.reg_cr,label,choiceInd,Thresh, sigThresh,savefluofigpath);
        end
        
        % outcome regression mask
        if ~exist('outcomeRegData','var') % if choice regression mask not computed
            label = 'outcome';
            outcomeInd = 7;
            outcomeRegData = getRegSelData(reg1.reg_cr,label,outcomeInd,Thresh, sigThresh,savefluofigpath);
        end
        % save the result
        
         % c(n+1) regression mask
        if ~exist('cn_1RegData','var') % if choice regression mask not computed
            label = 'upcoming choice';
            outcomeInd = 2;
            outcomeRegData = getRegSelData(reg1.reg_cr,label,outcomeInd,Thresh, sigThresh,savefluofigpath);
        end
        
        % r(n+1) regression mask
        if ~exist('rn_1RegData','var') % if choice regression mask not computed
            label = 'upcoming reward';
            rn_1Ind = 6;
            rn_1RegData = getRegSelData(reg1.reg_cr,label,rn_1Ind,Thresh, sigThresh,savefluofigpath);
        end
        
        % c(n-1) regression mask
        if ~exist('cn__1RegData','var') % if choice regression mask not computed
            label = 'previous choice';
            cn__1Ind = 4;
            cn__1RegData = getRegSelData(reg1.reg_cr,label,cn__1Ind,Thresh, sigThresh,savefluofigpath);
        end
        
        % r(n-1) regression mask
        if ~exist('rn__1RegData','var') % if choice regression mask not computed
            label = 'previous reward';
            rn__1Ind = 8;
            rn__1RegData = getRegSelData(reg1.reg_cr,label,rn__1Ind,Thresh, sigThresh,savefluofigpath);
        end
        
        % x(n) regression mask
        if ~exist('xnRegData','var') % if choice regression mask not computed
            label = 'interaction';
            xnInd = 11;
            xnRegData = getRegSelData(reg1.reg_cr,label,xnInd,Thresh, sigThresh,savefluofigpath);
        end
        
         if ~exist('xn__1RegData','var')
         label = 'previous interaction';
            xnInd = 12;
            xn__1RegData = getRegSelData(reg1.reg_cr,label,xnInd,Thresh, sigThresh,savefluofigpath);
         end
         
        % average reward regression mask
        if ~exist('ave_rRegData','var') % if choice regression mask not computed
            label = 'average reward';
            ave_rInd = 14;
            ave_rRegData = getRegSelData(reg1.reg_cr,label,ave_rInd,Thresh, sigThresh,savefluofigpath) ;
        end
        
        % cumulative reward regression mask
        if ~exist('cum_rRegData','var') % if choice regression mask not computed
            label = 'cumulative reward';
            cum_rInd = 15;
            cum_rRegData = getRegSelData(reg1.reg_cr,label,cum_rInd,Thresh, sigThresh,savefluofigpath) ;
        end
        
        %% regression 2-----------------------------------------------------------------------
        if ~exist('dQRegData','var') % if choice regression mask not computed
            label = 'delta Q';
            dQInd = 8;
            dQRegData = getRegSelData(reg2.reg_cr,label,dQInd,Thresh, sigThresh,savefluofigpath) ;
        end
        
        %if ~exist('chosenQRegData','var') % if choice regression mask not computed
            label = 'chosen Q';
            chosenQInd = 9;
            chosenQRegData = getRegSelData(reg2.reg_cr,label,chosenQInd,Thresh, sigThresh,savefluofigpath) ;
        %end
        
        if ~exist('dKRegData','var') % if choice regression mask not computed
            label = 'delta K';
            dKInd = 10;
            dKRegData = getRegSelData(reg2.reg_cr,label,dKInd,Thresh, sigThresh,savefluofigpath) ;
        end
        
        if ~exist('chosenKRegData','var') % if choice regression mask not computed
            label = 'chosen K';
            chosenKInd = 11;
            chosenKRegData = getRegSelData(reg2.reg_cr,label,chosenKInd,Thresh, sigThresh,savefluofigpath) ;
        end
        
        %% regression 3-----------------------------------------------------------------------
        
        if ~exist('RPERegData','var') % if choice regression mask not computed
            label = 'RPE';
            RPEInd = 6;
            RPERegData = getRegSelData(reg3.reg_cr,label,RPEInd,Thresh, sigThresh,savefluofigpath) ;
        end
        
        if ~exist('CKERegData','var') % if choice regression mask not computed
            label = 'CKE';
            CKEInd = 8;
            CKERegData = getRegSelData(reg3.reg_cr,label,CKEInd,Thresh, sigThresh,savefluofigpath) ;
        end
        
        %% save the results
        save(saveregmaskpath,'outcomeRegData','choiceRegData','cn_1RegData','cn__1RegData',...
                             'rn_1RegData','rn__1RegData', 'xnRegData', 'xn__1RegData', 'ave_rRegData', 'cum_rRegData',...
                             'dQRegData', 'chosenQRegData','dKRegData','chosenKRegData',...
                             'RPERegData', 'CKERegData')
    end
    
        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
    toc
end
end