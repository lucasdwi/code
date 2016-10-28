function [] = eventCollate(sdir,searchStr)
%% Get file structure for each searchStr
[fileStruct] = fileSearch(sdir,searchStr);
%eventLabel = [];
eventTimes = cell(numel(fileStruct));
timeMin = cell(numel(fileStruct));
timeMax = cell(numel(fileStruct));
for fsi = 1:numel(fileStruct)
    for fi = 1:size(fileStruct{fsi})
        load([sdir{1},fileStruct{fsi}(fi).name],'eventTs','LFPTs');
        indO = logicFind(0,cellfun('isempty',strfind(eventTs.label,'Orientation')),'==');
        indRs = logicFind(0,cellfun('isempty',strfind(eventTs.label,'Rest (Start)')),'==');
        indRe = logicFind(0,cellfun('isempty',strfind(eventTs.label,'Rest (End)')),'==');
        indAs = logicFind(0,cellfun('isempty',strfind(eventTs.label,'Approach (Start)')),'==');
        indAe = logicFind(0,cellfun('isempty',strfind(eventTs.label,'Approach (End)')),'==');
        indBs = logicFind(0,cellfun('isempty',strfind(eventTs.label,'Binge (Start)')),'==');
        indBe = logicFind(0,cellfun('isempty',strfind(eventTs.label,'Binge (End)')),'==');
        eventTimes{fsi}{fi,1} = eventTs.t{1,indO};
        eventTimes{fsi}{fi,2} = eventTs.t{1,indAs};
        eventTimes{fsi}{fi,3} = eventTs.t{1,indAe};
        eventTimes{fsi}{fi,4} = eventTs.t{1,indBs};
        eventTimes{fsi}{fi,5} = eventTs.t{1,indBe};
        eventTimes{fsi}{fi,6} = eventTs.t{1,indRs};
        eventTimes{fsi}{fi,7} = eventTs.t{1,indRe};
        timeMin{fsi}(fi) = LFPTs.tvec(1);
        timeMax{fsi}(fi) = LFPTs.tvec(end);
        %timeStep{fsi}(fi) = abs(0.0005-timeMin{fsi}(fi))/(timeMax{fsi}(fi)-timeMin{fsi}(fi));
        %tvec{fsi,fi} = LFPTs.tvec;
    end
end
%% Get percent of time at rest
for fsi = 1:numel(fileStruct)
    for fi = 1:size(fileStruct{fsi})
        restPerc(fsi,fi) = sum(eventTimes{fsi}{fi,7}-eventTimes{fsi}{fi,6})/(timeMax{fsi}(fi)-timeMin{fsi}(fi));
    end
end
%% Normalize and plot approach, binge, and rests
for fsi = 1:numel(eventTimes)
    for fi = 1:size(eventTimes{fsi},1)
        normApproach{fsi}{fi} = [((eventTimes{fsi}{fi,2}-timeMin{fsi}(fi))/(timeMax{fsi}(fi)-timeMin{fsi}(fi))),((eventTimes{fsi}{fi,3}-timeMin{fsi}(fi))/(timeMax{fsi}(fi)-timeMin{fsi}(fi)))];
        normBinge{fsi}{fi} = [((eventTimes{fsi}{fi,4}-timeMin{fsi}(fi))/(timeMax{fsi}(fi)-timeMin{fsi}(fi))),((eventTimes{fsi}{fi,5}-timeMin{fsi}(fi))/(timeMax{fsi}(fi)-timeMin{fsi}(fi)))];
        normRest{fsi}{fi} = [((eventTimes{fsi}{fi,6}-timeMin{fsi}(fi))/(timeMax{fsi}(fi)-timeMin{fsi}(fi))),((eventTimes{fsi}{fi,7}-timeMin{fsi}(fi))/(timeMax{fsi}(fi)-timeMin{fsi}(fi)))];
    end
end
%% Plot
figure
for fsi = 1:numel(fileStruct)
    subplot(2,2,fsi)
    q = 0;
    for fi = 1:numel(fileStruct{fsi})
        hold on;
        for j = 1:size(normApproach{fsi}{fi},1)
           h1 = plot([normApproach{fsi}{fi}(j,1) normApproach{fsi}{fi}(j,2)],[1+q 1+q],'-k');
        end
        for j = 1:size(normBinge{fsi}{fi},1)
           h2 = plot([normBinge{fsi}{fi}(j,1) normBinge{fsi}{fi}(j,2)],[1+q 1+q],'-r');
        end
        for j = 1:size(normRest{fsi}{fi},1)
            h3 = plot([normRest{fsi}{fi}(j,1) normRest{fsi}{fi}(j,2)],[1+q 1+q],'-bl');
        end
        q = q + 1;
    end
    xlabel('Normalized Time'); xlim([0 1]);
    ylabel('Animal'); ylim([0 numel(fileStruct{fsi})+1]);
    title(['Behavior Structure ', searchStr{fsi}])
end
legend([h1,h2,h3],'Approach','Binge','Rest')
%% Using timesteps get probability distributions for behaviors
% Create absoulte behavior vector
behTimeVect = (0:0.00001:1);
% Turn start and stop times into intervals
for fsi = 1:numel(fileStruct)
    behInts{fsi} = zeros(numel(fileStruct{fsi}),size(behTimeVect,2));
    for fi = 1:numel(fileStruct{fsi})      
        for ai = 1:size(normApproach{fsi}{fi},1)
            behInts{fsi}(fi,nearest_idx2(normApproach{fsi}{fi}(ai,1),behTimeVect):nearest_idx2(normApproach{fsi}{fi}(ai,2),behTimeVect)) = 1;
        end
        for bi = 1:size(normBinge{fsi}{fi},1)
            behInts{fsi}(fi,nearest_idx2(normBinge{fsi}{fi}(bi,1),behTimeVect):nearest_idx2(normBinge{fsi}{fi}(bi,2),behTimeVect)) = 2;
        end
        for ri = 1:size(normRest{fsi}{fi},1)
            behInts{fsi}(fi,nearest_idx2(normRest{fsi}{fi}(ri,1),behTimeVect):nearest_idx2(normRest{fsi}{fi}(ri,2),behTimeVect)) = 3;
        end
    end
end
%%
for ii = 1:numel(behInts)
    for j = 1:size(behInts{ii},2)
        appRatio(ii,j) = sum(behInts{ii}(:,j)==1)/size(behInts{ii},1);
        bingeRatio(ii,j) = sum(behInts{ii}(:,j)==2)/size(behInts{ii},1);
        restRatio(ii,j) = sum(behInts{ii}(:,j)==3)/size(behInts{ii},1);
    end
    cross = find((bingeRatio(ii,:)-restRatio(ii,:))>=0);
    figure
    hold on
    plot(behTimeVect,appRatio(ii,:),'-k')
    plot(behTimeVect,bingeRatio(ii,:),'-r')
    plot(behTimeVect,restRatio(ii,:),'-bl')
    %plot([behTimeVect(cross(end)) behTimeVect(cross(end))],[0 1],'-g')
    title(['Percent of ',searchStr{ii},' Animals in Behavior'])
    xlabel('Normalized Time')
    ylabel('Percent of Animals')
    legend('Approach','Binge','Rest')
end
%% Plot binge structure across conditions
figure
hold on
for ii = 1:4
    plot(behTimeVect,bingeRatio(ii,:))
end
legend('Base','24','48','Chow')
%% Get amount of time spent in binge and rest
for c = 1:numel(eventTimes)
    for ri = 1:size(eventTimes{c},1)
        % Approach
        timeSpent{c}(1,ri) = sum(eventTimes{c}{ri,3}-eventTimes{c}{ri,2})/timeMax{c}(ri);
        % Binge
        timeSpent{c}(2,ri) = sum(eventTimes{c}{ri,5}-eventTimes{c}{ri,4})/timeMax{c}(ri);
        % Rest
        timeSpent{c}(3,ri) = sum(eventTimes{c}{ri,7}-eventTimes{c}{ri,6})/timeMax{c}(ri);
    end
end
%% Get pop stats for binge and rest
for c = 1:numel(timeSpent)
    popStats(c,1) = mean(timeSpent{c}(2,:));
    popStats(c,2) = std(timeSpent{c}(2,:));
    popStats(c,3) = mean(timeSpent{c}(3,:));
    popStats(c,4) = std(timeSpent{c}(3,:));
end
%%
allRest = horzcat(timeSpent{1}(3,:),timeSpent{2}(3,:),timeSpent{3}(3,:),timeSpent{4}(3,:));
allBinge = horzcat(timeSpent{1}(2,:),timeSpent{2}(2,:),timeSpent{3}(2,:),timeSpent{4}(2,:));
allGroup = {'Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Base','Dep24','Dep24','Dep24','Dep24','Dep24','Dep24','Dep24','Dep24','Dep24','Dep24','Dep24','Dep24','Dep48','Dep48','Dep48','Dep48','Dep48','Dep48','Dep48','Dep48','Dep48','Chow','Chow','Chow','Chow','Chow','Chow','Chow','Chow'};
%% T-Tests
cmb = nchoosek(1:4,2);
for c = 1:size(cmb,1)
    [hrest{c},prest{c}] = ttest2(timeSpent{cmb(c,1)}(3,:),timeSpent{cmb(c,2)}(3,:));
    [hbinge{c},pbinge{c}] = ttest2(timeSpent{cmb(c,1)}(2,:),timeSpent{cmb(c,2)}(2,:));
end
%% Boxplots of Rest and Binge time
figure
b1 = boxplot(allRest,allGroup);
% hold on
% for ii = 1:length(prest)
%     if prest{ii} <= 0.05
%         plot(cmb(ii,1):cmb(ii,2),repmat(0.7,1,size(cmb(ii,1):cmb(ii,2),2)),'-k')
%     end
% end
title('Proportion of Time Spent Resting')
ylabel('Proportion of Time')
figure
boxplot(allBinge,allGroup)
title('Proportion of Time Spent Bingeing')
ylabel('Proportion of Time')
%%
restMean = [];
for ii = 1:2:length(test)
   restMean = vertcat(restMean,mean([test(ii),test(ii+1)]));
end