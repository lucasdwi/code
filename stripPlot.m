function stripPlot(data1,data2,p,x,main)
%% Creates scatter plot comparing two datasets and uses p-value to plot significance
% INPUTS:
% data1 = data array 1, format: observation X variable
% data2 = data array 2, format: observation X variable
% p = row vector of p-values, format: each column corresponds to
%   statistical test output from comparing that column in data1 to data2
% x = x-axis labels; format: row vector of strings
% main = title; format: string in cell
%%
figure
hold on
for c = 1:size(data1,2)
    if p(c) <= 0.05
        % Plot data1
        plot([c-0.25],data1(:,c),'.k')
        thism1 = mean(data1(:,c));
        % If comparing to zero, change color of mean square for above 
        % (blue) or below (red) 0
        if isempty(data2)
            if thism1 < 0
                plot([c-0.25],thism1,'sr')
            else if thism1 >0
                    plot([c-0.25],thism1,'sb')
                end
            end
        % If comparing two data sets, then just use red square for mean
        else
            plot([c-0.25],thism1,'sr')
        end
        % Plot data2
        if ~isempty(data2)
            plot([c+0.25],data2(:,c),'.','color',[.5 .5 .5])
            thism2 = mean(data2(:,c));
            plot([c+0.25],thism2,'sr')
        end
    end
end
% If comparing to zero, plot zero line
if isempty(data2)
   plot([1 size(data1,2)],[0 0],'k') 
end
% Using string vector of varaibles names for x-axis label
set(gca,'XTick',logicFind(0.05,p,'<='),'XTickLabels',x(logicFind(0.05,p,'<=')),'XTickLabelRotation',90)
% Set title
title(main)
