%% Box plots and T-Tests of just models with significant p-values
coreGroup = {'coreRespond','coreStrict'};
shellGroup = {'shellRespond','shellStrict'};
for i = 1:size(sigModels.Base,1)
    % if all and shell, or either shellrespond/strict
    if (strcmp(models.Base{2,sigModels.Base.ModelNumber(i)},'all') && strncmpi(models.Base{4,sigModels.Base.ModelNumber(i)},'s',1)) || strncmpi(models.Base{2,sigModels.Base.ModelNumber(i)},'s',1)
        thisGroup = shellGroup;
        % if all and core, or either corerespond/strict
    else if (strcmp(models.Base{2,sigModels.Base.ModelNumber(i)},'all') && strncmpi(models.Base{4,sigModels.Base.ModelNumber(i)},'c',1)) || strncmpi(models.Base{2,sigModels.Base.ModelNumber(i)},'c',1)
            thisGroup = coreGroup;
        end
    end        
    thisVar = cellstr(models.Base(3,sigModels.Base.ModelNumber(i)));
    % Setup data matrix
    data{i} = [T.Base.(thisVar{1})(T.Base.(thisGroup{1}) ~= 1);T.Base.(thisVar{1})(T.Base.(thisGroup{1}) == 1);T.Base.(thisVar{1})(T.Base.(thisGroup{2}) == 1)];
    len1 = length(T.Base.(thisVar{1})(T.Base.(thisGroup{1}) ~= 1));
    len2 = length(T.Base.(thisVar{1})(T.Base.(thisGroup{1}) == 1));
    len3 = length(T.Base.(thisVar{1})(T.Base.(thisGroup{2}) == 1));
    data{i}(1:len1,2) = 1; data{i}(len1+1:len1+len2,2) = 2; data{i}(len1+len2+1:len1+len2+len3,2) = 3;
    % T-tests
    [h{i,1},p{i,1}] = ttest2(data{i}(data{i}(:,2)==1,1),data{i}(data{i}(:,2)==2,1));       
    [h{i,2},p{i,2}] = ttest2(data{i}(data{i}(:,2)==1,1),data{i}(data{i}(:,2)==3,1));
    [h{i,3},p{i,3}] = ttest2(data{i}(data{i}(:,2)==2,1),data{i}(data{i}(:,2)==3,1));
end
%% Boxplots and T-tests of all models
coreGroup = {'coreRespond','coreStrict'};
shellGroup = {'shellRespond','shellStrict'};
for i = 1:size(models.Base,2)
    % if all and shell, or either shellrespond/strict
    if (strcmp(models.Base{2,i},'all') && strncmpi(models.Base{4,i},'s',1)) || strncmpi(models.Base{2,i},'s',1)
        thisGroup = shellGroup;
        % if all and core, or either corerespond/strict
    else if (strcmp(models.Base{2,i},'all') && strncmpi(models.Base{4,i},'c',1)) || strncmpi(models.Base{2,i},'c',1)
            thisGroup = coreGroup;
        end
    end        
    thisVar = cellstr(models.Base(3,i));
    % Setup data matrix
    data{i} = [T.Base.(thisVar{1})(T.Base.(thisGroup{1}) ~= 1);T.Base.(thisVar{1})(T.Base.(thisGroup{1}) == 1);T.Base.(thisVar{1})(T.Base.(thisGroup{2}) == 1)];
    len1 = length(T.Base.(thisVar{1})(T.Base.(thisGroup{1}) ~= 1));
    len2 = length(T.Base.(thisVar{1})(T.Base.(thisGroup{1}) == 1));
    len3 = length(T.Base.(thisVar{1})(T.Base.(thisGroup{2}) == 1));
    data{i}(1:len1,2) = 1; data{i}(len1+1:len1+len2,2) = 2; data{i}(len1+len2+1:len1+len2+len3,2) = 3;
    % T-tests
    [h{i,1},p{i,1}] = ttest2(data{i}(data{i}(:,2)==1,1),data{i}(data{i}(:,2)==2,1));       
    [h{i,2},p{i,2}] = ttest2(data{i}(data{i}(:,2)==1,1),data{i}(data{i}(:,2)==3,1));
    [h{i,3},p{i,3}] = ttest2(data{i}(data{i}(:,2)==2,1),data{i}(data{i}(:,2)==3,1));
end
%% Plots
% Check if any of comparisons reach significance, if so plot
clear b;
for i = 1:size(p,1)    
    if ((p{i,1}<=0.05) + (p{i,2}<=0.05) + (p{i,3}<=0.05)) >= 1
        figure;
        b{i} = boxplot(data{i}(:,1),data{i}(:,2),'Labels',{'Non-Responders','Responders','Strict'});
        yt = get(gca,'Ytick'); xt = get(gca,'Xtick');
        yTick = abs(yt(2)-yt(1));
        yMax = [get(b{i}(3,1),'YData');get(b{i}(3,2),'YData');get(b{i}(3,3),'YData')];
        hold on;
        for j = 1:3
            bar = {[1 2],[1 3],[2 3]};
            if p{i,j} <= 0.05
                if j == 1 || j == 3
                    plot(xt(bar{j}),[1 1]*max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+0.5*yTick,'-k',mean(xt(bar{j})),max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+0.8*yTick,'*k');
                    text(mean(xt(bar{j}))+0.1,max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+0.8*yTick,strcat('p \approx ',num2str(round(p{i,j},4))));
                end
                if j == 2
                    plot(xt(bar{j}),[1 1]*max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+yTick,'-k',mean(xt(bar{j})),max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+1.3*yTick,'*k');
                    text(mean(xt(bar{j}))+0.1,max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+1.3*yTick,strcat('p \approx ',num2str(round(p{i,j},4))));
                end
            end
        end
        lim = ylim;
        if lim(2) < max(yMax(:,1))+yTick*1.5
            ylim([lim(1) lim(2)+yTick*1.5]);
        end
        title(models.Base{4,i}); ylabel(models.Base(3,i)); 
        %title(models.Base{4,sigModels.Base.ModelNumber(i)}); ylabel(models.Base(3,sigModels.Base.ModelNumber(i)));
    end
 end