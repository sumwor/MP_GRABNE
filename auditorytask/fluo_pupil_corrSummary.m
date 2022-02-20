function fluo_pupil_corrSummary(dataIndex)

% summarize cross-correlation results of pupil and GRAB signal

nFiles = size(dataIndex,1);
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);

maxLagSum = [];
maxCorrSum = [];
for ii = 1:nFiles
     savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
     savefilename = fullfile(savematpath,'xcorr-pupil-fluo.mat');
     %fn_corr = dir(savefilename)
     if exist(savefilename)
         load(savefilename);
           % plot the cross-correlation with lags
         date = num2str(dataIndex.DateNumber(ii));
          savefigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_crosscorrelation']);
          if ~exist(savefigpath)
              mkdir(savefigpath);
          end
          
         maxCorr = zeros(1, length(c));
         maxLag = zeros(1, length(c));
         for rr = 1:length(c)
             % get the maximum correlation and lags
             corrInRange = c{rr}(newlag<=20 & newlag > -20);
             lagInRange = newlag(newlag<=20 & newlag > -20);
             [~,I] = max(abs(corrInRange));
             maxCorr(rr) = corrInRange(I);
             maxLag(rr) = lagInRange(I);
             
             aCorrF = autoCorrFluo{rr}(newlagFluo<=100 & newlagFluo>=0);
             aLagF= newlagFluo(newlagFluo<=100 & newlagFluo>=0);
             
             
             % plot cross- and auto-correlation
             %if ~exist(fullfile(savefigpath,['CrossCorrelationRoi-',num2str(length(c)),'.fig']))
                figure;stem(lagInRange,corrInRange);
                print(gcf,'-dpng',fullfile(savefigpath,['CrossCorrelationRoi-',num2str(rr)]));
                saveas(gcf, fullfile(savefigpath,['CrossCorrelationRoi-',num2str(rr)]), 'fig');
                
                figure;stem(aLagF,aCorrF);
                print(gcf,'-dpng',fullfile(savefigpath,['autoCorrelationRoi-',num2str(rr)]));
                saveas(gcf, fullfile(savefigpath,['autoCorrelationRoi-',num2str(rr)]), 'fig');
                close all;
             %end
         end
         aCorrP = autoCorrPup(newlagPup<=100 & newlagPup>=0);
         aLagP = newlagPup(newlagPup<=100 & newlagPup>=0);
            figure;stem(aLagP,aCorrP);
                print(gcf,'-dpng',fullfile(savefigpath,['autoCorrelationPup']));
                saveas(gcf, fullfile(savefigpath,['autoCorrelationPup']), 'fig');
                close all;
        
         % plot the maxcorr and maxlag in 2d
         tlabel = 'Correlation coefficient';
         colorRange = [-0.3 0.3];
         plot_xcorr(maxCorr,tlabel,colorRange,savefigpath);
         % max lag
          tlabel = 'lags';
         colorRange = [-20 20];
         date = num2str(dataIndex.DateNumber(ii));
         savefigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
         plot_xcorr(maxLag,tlabel,colorRange,savefigpath)
         % plot distributions of correlation and lag within single session
         
         
         close all;
         maxCorrSum = cat(1,maxCorrSum, maxCorr);
         maxLagSum= cat(1, maxLagSum, maxLag);
    
   
     
     %% summarize the lags and xcorr
    
     end
end

 display('Sum');
end

