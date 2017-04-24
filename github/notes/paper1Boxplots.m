% c = 5 % Binge size food dep
% c = 6; % CPP
c = 8; % Rest ratio
figure;
plot([1],table2array(T.Base(nonCI,c)).*100,'ko')
hold on; plot([2],table2array(T.Base(respCI,c)).*100,'ko')
hold on; plot([3],table2array(T.Base(strictCI,c)).*100,'ko')
hold on; plot([.75 1.25],[nanmean(table2array(T.Base(nonCI,c))) nanmean(table2array(T.Base(nonCI,c)))].*100,'-r')
hold on; plot([1.75 2.25],[nanmean(table2array(T.Base(respCI,c))) nanmean(table2array(T.Base(respCI,c)))].*100,'-r')
hold on; plot([2.75 3.25],[nanmean(table2array(T.Base(strictCI,c))) nanmean(table2array(T.Base(strictCI,c)))].*100,'-r')
% title('Core DBS Outcome vs. Food Dep. Binge Size Change')
% title('Shell DBS Outcome vs. CPP')
title('Core DBS Outcome vs. Rest Ratio')
% ylabel('% Change in Binge Size')
% ylabel('% Change in NPC Time')
ylabel('% of Time at Rest')
xtick([1,2,3]); set(gca,'XTickLabel',{'Non-Responder','All','Strict'}); xlim([.5 3.5])
%%
%hold on; plot([1 3], [64.5 64.5],'-k')
%hold on; plot([2], [65.5],'k*')
%hold on; text(2.1,65.5,['p \approx ',num2str(0.0027)])
hold on; plot([1 2], [64.5 64.5],'-k')
hold on; plot([1.5], [65.5],'k*')
hold on; text(1.6,65.5,['p \approx ',num2str(0.0112)])