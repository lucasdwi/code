function [p,pAdj] = anovaBox(x,group,varName,varUnit,plotType)
%% Runs a 1-way ANOVA test and plots boxplots with significance.
% INPUTS:
% x = data; format: column vector
% group = group assignment matching data column vector; format: vector of
%   strings or numerals
% varName = variable being tested, used for title; format: string
% varUnit = unit of variable being tested, used for y-axis label; format:
%   string
% plotType = type of plot to use; format: 'box' or 'dot'

% OUTPUTS:
% p = group difference p-value from anova1; format: numeric
% pAdj = multiple comparison corrected (Bonferroni) p-values; format: n x
%   1 column vector of p-values, where n = number of comparisons
%   (nchoosek(number of groups,2))

% * = p<0.05
% ** = p<0.01
% *** = p<0.001
%% Written by LLD 2017; Uses sigstar by Rob Campbell - CSHL 2013
%%
% Run ANOVA
% If group not supplied, run with each column as group
if isempty(group)
    [p,~,stats] = anova1(x,[],'off');
else
    [p,~,stats] = anova1(x,group,'off');
end
% Run multiple corrections (Bonferroni) on all pairs of means
[mctbl] = multcompare(stats,'CType','bonferroni','display','off');
% Combine group labels for sigstar
sigGroups = cell(nchoosek(numel(unique(group)),2));
for ii = 1:size(mctbl,1)
    sigGroups{ii} = [mctbl(ii,1),mctbl(ii,2)];
end
% Set-up plot
if strcmpi(plotType,'box')
    figure;
    boxplot(x,group)
    hold on
    for ii = 1:numel(unique(group))
        plot(ii,stats.means(ii),'rs')
    end
elseif strcmpi(plotType,'dot')
    [xSpace] = stripPlot(x,group);
end
% Plot means as red squares

title([varName,' across Groups'])
ylabel([varName,' ',varUnit]);
% Find indices of mean comparisons with significant differences
sigInds = logicFind(0.05,mctbl(:,6),'<=');
% Plot significance bars with sigstar
if strcmpi(plotType,'box')
    sigstar(sigGroups(sigInds),mctbl(sigInds,6))
elseif strcmpi(plotType,'dot')
    %adjSigGroup = cellfun(@minus,sigGroups(sigInds),{1},'UniformOutput',0);
    sigstar(xSpace(cell2mat(sigGroups(sigInds))),mctbl(sigInds,6))
end
pAdj = mctbl(:,6);