function MP_GRAB_temporalCorrSummary(dataIndex, savematsumpath)

% calculate temporal correlation of different varible coefficients within
% grids

% check if coefficient sign changes over time
% for eachgrid, calculate the number of + time point and - time point
% then average

nFiles = size(dataIndex,1);


%% get clustering results, average similar clusters, also rerun clustering on all grids together
for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    %load(fullfile(fn_beh.folder,fn_beh.name));


    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1

        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end

        cd(savefluofigpath);


        %load(fullfile(fn_fluo.folder,fn_fluo.name));

        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        
        % outcome
        tlabel1 = 'Outcome coefficient';
        saveregpath =  fullfile(savematpath,[tlabel1,' cluster.mat']);  % file path to save the results
        if exist(saveregpath)
            load(saveregpath)
        end

    end
end

%% get lag
% cnLagCorr = zeros(1,nFiles); cnLagP =  zeros(1,nFiles); cnDis = []; cnCorr = []; cnLag = [];
% rnLagCorr = zeros(1,nFiles); rnLagP = zeros(1,nFiles); rnDis = []; rnCorr = []; rnLag = [];
% xnLagCorr = zeros(1,nFiles); xnLagP = zeros(1,nFiles); xnDis = []; xnCorr = []; xnLag = [];
% 
% 
% for ii = 1:nFiles
% 
%     % load behavior files
%     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
%     %load(fullfile(fn_beh.folder,fn_beh.name));
% 
% 
%     % load fluorescent files
%     fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
%     if length(fn_fluo) == 1
% 
%         savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
%         if ~exist(savefluofigpath,'dir')
%             mkdir(savefluofigpath);
%         end
% 
%         cd(savefluofigpath);
% 
% 
%         %load(fullfile(fn_fluo.folder,fn_fluo.name));
% 
%         % make folders to save analysis and plots
%         savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
%         if ~exist(savematpath,'dir')
%             mkdir(savematpath);
%         end
% 
%         saveregpath = fullfile(savematpath,'regressionVarTempCorr.mat');  % file path to save the results
%         if exist(saveregpath)
%             load(saveregpath)
%         end
% 
%         % get the data
%         % scatter plot of distance from reference grid and lag
%         
%         % get distance
%         cnDistance = []; cnlag = []; cnCoeff = [];
%         for tt = 1:size(choicetempData.tempCorrLag,1)
%            for uu = 1:size(choicetempData.tempCorrLag,2)
%                if ~isnan(choicetempData.tempCorrLag(tt,uu))
%                    cnDistance = [cnDistance,sqrt((tt-choicetempData.maxX)^2+(uu-choicetempData.maxY)^2)];
%                    cnlag = [cnlag, choicetempData.tempCorrLag(tt,uu)];
%                    cnCoeff = [cnCoeff, choicetempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%         
%         [R,P] = corrcoef(cnDistance,abs(cnlag));
%         if ~isnan(R)
%             cnLagCorr(ii) = R(1,2); cnLagP(ii) = P(1,2); cnDis = [cnDis,cnDistance]; cnCorr = [cnCorr, cnCoeff]; cnLag = [cnLag, cnlag];
%         end
%             %figure;scatter(cnDistance,cnLag);
%         
%         rnDistance = []; rnlag = []; rnCoeff = [];
%         for tt = 1:size(outcometempData.tempCorrLag,1)
%            for uu = 1:size(outcometempData.tempCorrLag,2)
%                if ~isnan(outcometempData.tempCorrLag(tt,uu))
%                    rnDistance = [rnDistance,sqrt((tt-outcometempData.maxX)^2+(uu-outcometempData.maxY)^2)];
%                    rnlag = [rnlag, outcometempData.tempCorrLag(tt,uu)];
%                    rnCoeff = [rnCoeff, outcometempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%            [R,P] = corrcoef(rnDistance,abs(rnlag));
%             if ~isnan(R)
%         rnLagCorr(ii) = R(1,2); rnLagP(ii) = P(1,2); rnDis = [rnDis,rnDistance]; rnCorr = [rnCorr, rnCoeff]; rnLag = [rnLag, rnlag];
%             end
% 
%         xnDistance = []; xnlag = []; xnCoeff = [];
%         for tt = 1:size(xntempData.tempCorrLag,1)
%            for uu = 1:size(xntempData.tempCorrLag,2)
%                if ~isnan(xntempData.tempCorrLag(tt,uu))
%                    xnDistance = [xnDistance,sqrt((tt-xntempData.maxX)^2+(uu-xntempData.maxY)^2)];
%                    xnlag = [xnlag, xntempData.tempCorrLag(tt,uu)];
%                    xnCoeff = [xnCoeff, xntempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%            [R,P] = corrcoef(xnDistance,abs(xnlag));
% 
%             if ~isnan(R)
%         xnLagCorr(ii) = R(1,2); xnLagP(ii) = P(1,2); xnDis = [xnDis,xnDistance]; xnCorr = [xnCorr, xnCoeff]; xnLag = [xnLag, xnlag];
%             end
%     end
% end    
% 
% % save the results
% savematdir = fullfile(savematsumpath, 'tempCorr');
% if ~exist(savematdir)
%     mkdir(savematdir)
% end
% cnTemp.dis = cnDis; cnTemp.corr = cnCorr; cnTemp.lag = cnLag; cnTemp.corrp = cnLagP; cnTemp.corrCoef = cnLagCorr;
% rnTemp.dis = rnDis; rnTemp.corr = rnCorr; rnTemp.lag = rnLag; rnTemp.corrp = rnLagP; rnTemp.corrCoef = rnLagCorr;
% xnTemp.dis = xnDis; xnTemp.corr = xnCorr; xnTemp.lag = xnLag; xnTemp.corrp = xnLagP; xnTemp.corrCoef = xnLagCorr;
% save(fullfile(savematdir,'tempCorrLag.mat'),'cnTemp','rnTemp','xnTemp');
% 
% cnLagCorrSig = cnLagCorr(cnLagP<0.05);
% rnLagCorrSig = rnLagCorr(rnLagP<0.05);
% xnLagCorrSig = xnLagCorr(xnLagP<0.05);
% 
%  [p,h,stats] = signrank(cnLagCorrSig,0,'tail','right')
%   [p,h,stats] = signrank(rnLagCorrSig,0,'tail','right')
%    [p,h,stats] = signrank(xnLagCorrSig,0,'tail','right')
% 
% figure;scatter(cnDis, abs(cnLag),'.');
% [R,P] = corrcoef(cnDis,abs(cnLag))
% 
% figure;scatter(rnDis, abs(rnLag),'.');
% [R,P] = corrcoef(rnDis,abs(rnLag))
% 
% figure;scatter(xnDis, abs(xnLag),'.');
% [R,P] = corrcoef(xnDis,abs(xnLag))