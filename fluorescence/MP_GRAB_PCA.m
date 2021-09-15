function MP_GRAB_PCA(dataIndex)

%% pca analysis
% dPCA
%  also try t-SNE?

% separate the trials? or not?

nFiles = size(dataIndex,1);

%% go through every session, run analysis
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    load(fullfile(fn_fluo.folder,fn_fluo.name));
    
    fn_ROI = dir(fullfile(dataIndex.LogFilePath{ii},'cell*.mat'));
    
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
    fn_PCA = dir(fullfile(savematpath,'pca.mat'));  % regression for fluo change
    
    if isempty(fn_PCA)
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
        %% get meam intensity
        meanProjPath = dir(fullfile(dataIndex.BehPath{ii},'*Mean*.tif'));
        meanProj = loadtiffseq([],fullfile(meanProjPath.folder,meanProjPath.name));
        
        % get the mean intensity in ROI
        meanIntensity = zeros(1, size(fn_ROI,1));
        %         varIntensity = zeros(1,size(fn_ROI,1));
        for uu = 1:size(fn_ROI,1)
            roi = load(fullfile(fn_ROI(uu).folder,fn_ROI(uu).name));
            meanIntensity(uu) = mean(meanProj(roi.bw==1));
            %             varIntensity(uu) = nanstd(cells.cellf{uu});
        end
        
        % save figure path
        %% run PCA on dFF of all 196 ROIs
        % get the ROI data
        %         pca_observe = zeros(length(cells.dFF{1}),length(cells.dFF));
        %         for jj = 1:length(cells.dFF)
        %             pca_observe(:,jj) = cells.dFF{jj};
        %         end
        %         [coeff,score,latent,tsquared,explained,mu] = pca(pca_observe);
        %
        %         figure;
        % %         biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4'});
        % %             close all;
        %         scatter3(score(:,1),score(:,2),score(:,3))
        %         axis equal
        %         xlabel('1st Principal Component')
        %         ylabel('2nd Principal Component')
        %         zlabel('3rd Principal Component')
        
        %% dPCA
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = [];  %confidence interval
        params.minNumTrial = 50;
        
        for uu=1:numel(cells.dFF)
            %numROI = sum(meanIntensity(uu) > prctile(meanIntensity,45));
            %if meanIntensity(uu) > prctile(meanIntensity,45)
                psth_output=[];
                
                fieldname{1}={'reward'};
                fieldname{2}={'noreward'};
%                 fieldname{1}={'left','reward'};
%                 fieldname{2}={'left','noreward'};
                
%                 fieldname{3}={'right','reward'};
%                 fieldname{4}={'right','noreward'};
                
                for kk=1:numel(fieldname)
                    trialMask = getMask(trials,fieldname{kk});
                    out = get_psth( cells.dFF{uu}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
                    if kk == 1
                        psth_output = zeros(length(cells.cellf),length(out.t),numel(fieldname));
                    end
                    psth_output(uu,:,kk) = out.signal;
                end
                
                
            %end
        end
        
        % running dPCA
        comps = 20;
        %maxstep = 1e5;
        
        W = dpca(psth_output,comps,[],[]);
    end
end
end


