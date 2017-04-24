%% Open behavior data
load('C:\Users\Lucas\Desktop\GreenLab\data\behaviorData.mat')
cd('C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed')
%% Open all baseline files and store individual and average binge size
files = []; files = dir('*Base*');
for f = 1:size(files,1)
    load(files(f).name);
    fRow = find(strcmp(files(f).name(1:end-4),table2cell(behaviorData(:,1)))==1);
    bingeSize = table2array(behaviorData(fRow,4));
    %avgBinge = table2array(behaviorData(fRow,5));
    save(files(f).name,'pl2','lfpchan','LFPTs','adfreq','ts','fn','ad','n','WBchan','TimeSampEr','eventTs','bingeSize');
    clearvars -except behaviorData files f
end
% Open all FoodDep (24 and 48) and RegChow and store % change from baseline and bingeSize
files = []; files = vertcat(dir('*Dep*'),dir('*Chow*'));
for f = 1:size(files,1)
    load(files(f).name);
    fRow = find(strcmp(files(f).name(1:end-4),table2cell(behaviorData(:,1)))==1);
    bingeSize = table2array(behaviorData(fRow,4));
    %bingeChangeAvg = table2array(behaviorData(fRow,6));
    %bingeChangeOne = table2array(behaviorData(fRow,7));
    save(files(f).name,'pl2','lfpchan','LFPTs','adfreq','ts','fn','ad','TimeSampEr','eventTs','n','WBchan','bingeSize');
    clearvars -except behaviorData files f
end 