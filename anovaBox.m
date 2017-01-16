function [p] = anovaBox(x,group,varName,varUnit)
% Run ANOVA
[p,tbl,stats] = anova1(x,group,'off');
% Run multiple corrections (Bonferroni) on all pairs of means
[mctbl] = multcompare(stats,'CType','bonferroni','display','off');
% Combine group labels for sigstar
for ii = 1:size(mctbl,1)
    sigGroups{ii} = [mctbl(ii,1),mctbl(ii,2)];
end
% Set-up boxplot
figure; 
boxplot(x,group)
hold on
% Plot means as red squares
for ii = 1:length(group)
   plot(ii,stats.means(ii),'rs') 
end
title([varName,' across Conditions'])
ylabel([varName,' ',varUnit]);
% Find indices of mean comparisons with significant differences
sigInds = logicFind(0.05,mctbl(:,6),'<=');
% Plot significance bars with sigstar
sigstar(sigGroups(sigInds),mctbl(sigInds,6))