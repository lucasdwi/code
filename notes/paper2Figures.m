%% Run ANOVA and post hoc Bonferroni-corrected mean comparisons; plot boxplots with sig bars
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
% Run one way ANOVA
[p,tbl,stats] = anova1(bingeSizes,group,'off');
% Run multiple corrections (Bonferroni) on all pairs of means
[mctbl] = multcompare(stats,'CType','bonferroni','display','off');
% Combine group labels for sigstar
for ii = 1:size(mctbl,1)
    sigGroups{ii} = [mctbl(ii,1),mctbl(ii,2)];
end
% Set-up boxplot
figure; 
boxplot(bingeSizes,group)
hold on
% Plot means as red squares
for ii = 1:length(group)
   plot(ii,stats.means(ii),'rs') 
end
title('Binge Size across Conditions')
ylabel('Binge Size (gm)');
% Find indices of mean comparisons with significant differences
sigInds = logicFind(0.05,mctbl(:,6),'<=');
% Plot significance bars with sigstar
sigstar(sigGroups(sigInds),mctbl(sigInds,6))
%% Try predicting binge vs rest
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\H10BaseOct15_binge_vs_rest.mat')
%%
binge = [reshape(relPower.event1,1,48),reshape(coh{1,1}.rel,1,36];
rest = [reshape(relPower.event2,1,48),];