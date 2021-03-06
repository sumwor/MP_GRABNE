function plot_switch_hrside(input,tlabel)
% % plot_switch_hrside %
%PURPOSE:   Plot choice behavior around switches, when the high reward
%           side has changed
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   input:  Structure generated by choice_switch_hrside().
%   tlabel: Text string that will be put on top of the figure

%%
if numel(input)>1  %more than 1 data set, plot mean+-sem
    n=input{1}.n;
    for j=1:numel(input)
        temp_probh(:,j)=input{j}.probh;
        temp_probl(:,j)=input{j}.probl;
    end
    probh=nanmean(temp_probh,2);
    probl=nanmean(temp_probl,2);
    probh_sem=nanstd(temp_probh,[],2)./sqrt(numel(input));
    probl_sem=nanstd(temp_probl,[],2)./sqrt(numel(input));
else        %plot the 1 data set
    n=input.n;
    probh=input.probh;
    probl=input.probl;
end

figure;
subplot(1,2,1); hold on;
plot(n,probh,'.-','MarkerSize',30,'Linewidth',3,'Color','k');
plot(n,probl,'.-','MarkerSize',30,'Linewidth',3,'Color','m');
plot([0 0],[0 1],'k--','LineWidth',2);
if numel(input)>1
    for l=1:numel(n)
        plot(n(l)*[1 1],probh(l)+probh_sem(l)*[-1 1],'-','LineWidth',3,'Color','k');
        plot(n(l)*[1 1],probl(l)+probl_sem(l)*[-1 1],'-','LineWidth',3,'Color','m');
    end
end
ylabel('P_{choice}');
xlabel('Trials from block switch');
legend('High reward->low reward','Low reward->high reward');
axis([n(1) n(end) 0 1]);
title(tlabel);

print(gcf,'-dpng','switches_hrside');    %png format
saveas(gcf, 'switches_hrside', 'fig');

end


    