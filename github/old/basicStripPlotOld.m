function basicStripPlot(data,group,meanMark)
%%
% either plot every column or plot group divided columns

%% Compare each column across groups
% NOTE: to generalize beyond binary group, use unique on group, then cycle
% through those 
if ~isempty(group)
    for di = 1:size(data,2)
        figure
        hold on
        % Plot group 1
        plot(1,data(group==1,di),'.k','MarkerSize',10)
        % Plot group 1 mean
        if strcmpi(meanMark,'s')
            plot(1,nanmean(data(group==1,di)),'rs','MarkerFaceColor','r')
        else if strcmpi(meanMark,'l')
                plot([0.75 1.25],[nanmean(data(group==1,di)) nanmean(data(group==1,di))],'-k','LineWidth',2,'Color',[0.5 0.5 0.5])
            end
        end
        % Plot group 2
        plot(2,data(group==0,di),'.k','MarkerSize',10)
        % Plot group 2 mean
        if strcmpi(meanMark,'s')
            plot(2,nanmean(data(group==0,di)),'rs','MarkerFaceColor','r')
        else if strcmpi(meanMark,'l')
                plot([1.75 2.25],[nanmean(data(group==0,di)) nanmean(data(group==0,di))],'-k','LineWidth',2,'Color',[0.5 0.5 0.5])
            end
        end
        xlim([0 3])
    end
end
%% Plot all columns together
if isempty(group)
    figure
    hold on
    for di = 1:size(data,2)
        plot(di,data(:,di),'.k','MarkerSize',10)
        if strcmpi(meanMark,'s')
            plot(di,nanmean(data(group==1,di)),'rs','MarkerFaceColor','r')
        else if strcmpi(meanMark,'l')
                plot([di-0.25 di+0.25],[nanmean(data(group==1,di)) nanmean(data(group==1,di))],'-k','LineWidth',2,'Color',[0.5 0.5 0.5])
            end
        end
    end
end
