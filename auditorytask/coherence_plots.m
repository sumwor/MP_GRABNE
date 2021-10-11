function coherence_plots(dataIndex)
% plot spontaneous pupil & fluorescent dynamics

nFiles = size(dataIndex,1);
cxyAll_NE.coeff = [];
cxyAll_Ach.coeff = [];
for ii = 1:nFiles
    
    % load behavior files
     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii},['*',date(1:6),'*_pup.mat']);
    dff_name = fullfile(fn_beh.folder,'dff.mat');
    
    fn_pup = dir(pup_name);
    fn_dff = dir(dff_name);
   
    
    if length(fn_pup) == 1 & length(fn_dff) == 1
         savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        % coherence analysis
       
       
        
        savematname='coherence.mat';
        cohfilename = fullfile(savematpath,savematname);
        load(cohfilename);
        if ii <= 3
            cxyAll_NE.coeff = [cxyAll_NE.coeff; cxy.coeff];
        else
            cxyAll_Ach.coeff = [cxyAll_Ach.coeff; cxy.coeff];
        end
    end
        
end

% plot total of 588 grids
cxyAll_NE = getBootstrp(cxyAll_NE, 0, 0.05);
cxyAll_Ach = getBootstrp(cxyAll_Ach, 0, 0.05);        

c_NE = [0.8500 0.3250 0.0980];
c_Ach = [0.9290 0.6940 0.1250];
figure;



%patch([fc fliplr(fc)], [cxy.boothigh  fliplr(cxy.bootlow)], [0.7 0.7 0.7])

errorshade(fc(2:end),cxyAll_NE.bootlow(2:end),cxyAll_NE.boothigh(2:end),c_NE,0.2);
hold on;
plot(fc(2:end),cxyAll_NE.coeff_bootave(2:end),'-','Color',c_NE);

hold on;
errorshade(fc(2:end),cxyAll_Ach.bootlow(2:end),cxyAll_Ach.boothigh(2:end),c_Ach,0.2);
hold on;
plot(fc(2:end),cxyAll_Ach.coeff_bootave(2:end),'-','Color',c_Ach);

xlim([0.001, 1]);
set(gca,'box','off')
set(gca, 'XScale', 'log')

xlabel('Coherence (with pupil)');

print(gcf,'-dpng',['spon_pupil_NE_Ach']);
saveas(gcf, 'spon_pupil_NE_Ach', 'fig');
saveas(gcf, 'spon_pupil_NE_Ach', 'svg');
end