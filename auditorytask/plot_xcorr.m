function plot_xcorr(input,tlabel,colorRange,savefigpath)
% % plot_selectivity %
%PURPOSE:   Plot selectivity, based on PSTHs from two conditions
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   input:       Structure generated by calc_selectivity()
%   sortParam:
%       if #values = 2, e.g. [1 3], then sort the cells based on selectivity value in this time period
%       if #values = number of cells, e.g. [1 3 4 2 5], then sort the cells directly based on this ordering
%   tlabel:       Text to put as title of the plot.
%   xtitle:       Text to put as the label for x-axis.
%   colorRange:   Max and min values that define the range of the color scale
%
%OUTPUT ARGUMENTS
%   cellOrder:    Sorting order of the cells
%                 (could be fed back as input 'sortParam' for next call of
%                 this function)

%% Plotting parameters
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);

%% Setup
nCells=numel(input);   
    
%%
% plot the selectivity in 2-D at differnet time point
edgelength = sqrt(nCells);

corr = zeros(edgelength);

for cc = 1:nCells
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

    corr(Ind1, Ind2) = input(cc);

end
figure;
subplot(1,2,1)
image(corr,'CDataMapping','scaled');
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);      %normalize dF/F heatmap to max of all conditions
ylabel('Cells');
%xlabel(xtitle);


%make a color scale bar
subplot(3,20,55);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
title(tlabel);
set(gca,'YDir','normal');
set(gca,'XTick',[]);
print(gcf,'-dpng',fullfile(savefigpath,tlabel));
            saveas(gcf, fullfile(savefigpath,tlabel), 'fig');

          
       

end