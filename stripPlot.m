function [xSpace] = stripPlot(data,group,varNames)
%% Plots data by group using dots for each datapoint and squares for group 
% means.  

% INPUTS:
% data = data to be visualized; format: either column vector or matrix with
%   observation X variable
% group = group assignment per observation; format: column vector of equal
%   length as number of rows in data
% varNames = optional input for matrix data, names corresponding to columns
%   of matrix; format: string array
%%
% If vector, plots with comparisons across groups. 
% If matrix, two options: (1) all same group and plots each column as own
% distribution (2) separate groups - as indicated by group - and plots
% groups side by side for each column.

% Transpose data and group if given as row vectors
if isrow(data)
    data = data';
end
if isrow(group)
    group = group';
end
% Check that data dimensions and group dimensions match
if size(data,1) ~= length(group)
   error('Data and groups do match dimensionality! Number of data points does not match group labels.') 
end
% If groups are cells of strings, convert to numbers for plotting, but keep
% string version
if iscellstr(group)
    groupStr = group;
    uGroup = unique(group);
    group = zeros(size(groupStr,1),1);
    for iG = 1:numel(uGroup)
        gInds = logicFind(uGroup(iG),groupStr,'==');
        group(gInds,1) = iG;
    end
end
% Grab unique groups
uGroup = unique(group);
% Determine spacing for plotting
spacing = fliplr(1./(1:numel(uGroup)));
% Grab midpoint for labels
mid = (spacing(end)-spacing(1))/2;
% Set up figure
figure
for iC = 0:size(data,2)-1
    for iG = 1:numel(uGroup)
        % Plot group data
        plot(repmat(iC+spacing(iG),sum(group==uGroup(iG)),1),data(group==uGroup(iG),iC+1),'.k')
        hold on
        % Plot mean of group
        plot(iC+spacing(iG),mean(data(group==uGroup(iG),iC+1),'omitNaN'),'rs')
    end
end
xlim([0 spacing(end)+iC+0.5])
% % Use group labels, either string or numeral, as x labels
% xSpace = spacing(1):spacing(1):spacing(end)*(iC+1);
% if exist('groupStr','var')
%     set(gca,'XTick',xSpace,'XTickLabel',unique(groupStr))
% else
%     set(gca,'XTick',xSpace,'XTickLabel',uGroup)
% end
% % If data is matrix, add second x-axis above to delineate variables
% if exist('varNames','var') && size(data,2)
%    ax1 = gca;
%    axes('Position',ax1.Position,'XAxisLocation','top','Color','none','YTickLabel','','XLim',ax1.XLim,'XTick',1-mid:spacing(end)*(iC+1)-mid,'XTickLabel',varNames);
% end
