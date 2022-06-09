function MP_GRAB_tempComp(save_path_mat_NE, save_path_mat_ACh)


savematNE = fullfile(save_path_mat_NE, 'tempCorr','tempCorrLag.mat');
savematACh = fullfile(save_path_mat_ACh, 'tempCorr','tempCorrLag.mat');

NETemp = load(savematNE);
AChTemp = load(savematACh);

% ,'ViolinColor',colors,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]
choice.NE = NETemp.cnTemp.corrCoef;
choice.ACh = AChTemp.cnTemp.corrCoef;
figure;violinplot(choice);

outcome.NE = NETemp.rnTemp.corrCoef;
outcome.ACh = AChTemp.rnTemp.corrCoef;
figure;violinplot(outcome);

interaction.NE = NETemp.xnTemp.corrCoef;
interaction.ACh = AChTemp.xnTemp.corrCoef;
figure;violinplot(interaction);


[p,h]=ranksum(NETemp.cnTemp.corrCoef,AChTemp.cnTemp.corrCoef)
[p,h]=ranksum(NETemp.rnTemp.corrCoef,AChTemp.rnTemp.corrCoef)
[p,h]=ranksum(NETemp.xnTemp.corrCoef,AChTemp.xnTemp.corrCoef)