function plotMaskBoundary(mergeMask,Mask1,Mask2, label,savefluofigpath,isSig)

%% isSig: 
%       0:coefficient sign (1 for positive,-1 for negative)
%       1:pval
%       2:coefficient value 
if isSig==1  % signficance mask
    colorRange = [0 1];
    colors=cbrewer('seq','Reds',256);
    titleText = ['Average ', label, ' significance mask'];
elseif isSig == 0   
    colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    colorRange = [-1 1];
     titleText = ['Average ', label, ' coefficient sign mask'];
else   
    colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    bound = max(abs(max(max(nanmean(mergeMask,3)))),abs(min(min(nanmean(mergeMask,3)))));
    
    colorRange = [-abs(bound), abs(bound)];
     titleText = ['Average ', label, ' coefficient value mask'];
end


figure;
subplot(1,2,1)
image(nanmean(mergeMask,3),'CDataMapping','scaled');
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
title(titleText);
hold on;
%get a function to plot different Masks
% plotbound(Mask1);plotbound(Mask2);%plotbound(Mask3);%plotbound(Mask4);
% title([label,'-selecvitity-mask']);
subplot(3,20,60);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
print(gcf,'-dpng',fullfile(savefluofigpath,titleText));
saveas(gcf, fullfile(savefluofigpath,titleText), 'fig');
saveas(gcf, fullfile(savefluofigpath,titleText), 'svg');

close;
