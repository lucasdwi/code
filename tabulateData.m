function [T,allPredict,allResponse,allGroups] = tabulateData
files.Base = {'H10BaseSep27','H13BaseSep27','H14BaseSep29','H15BaseSep29','I1BaseNov9','I2BaseNov9','I3BaseNov9','I4BaseSep24','I6BaseSep18','I8BaseSep23','I11BaseOct30','I12BaseNov12'};
files.Dep = {'H10FoodDepSep25','H13FoodDepSep25','H14FoodDepSep28','H15FoodDepSep28','I1FoodDepNov13','I2FoodDep24Dec16','I3FoodDep24Nov2','I4FoodDepSep30','I6FoodDepSep30','I8FoodDepSep30','I11FoodDep24Dec16','I12FoodDepNov13'};
%% Extract data from all processed files
cond = fieldnames(files); %{'Base','Dep'};
PSDs = {};
Cohs = {};
ntrials = {};
for c = 1:length(cond)
    events = {'event1','event2'};
    for i = 1:length(files.(cond{c}))
        load(strcat('C:\Users\Lucas\Desktop\GreenLab\data\processed\',files.(cond{c}){i},'_processed.mat'));
        for j = 1:length(events)
            PSDs.(cond{c}).Overall{i,j} = psdTrls.(events{j}).Overall;
            PSDs.(cond{c}).Bands{i,j} = psdTrls.(events{j}).Avg;
            ntrials.(cond{c}){i,j} = length(psdTrls.(events{j}).Pow);
            Cohs.rel.(cond{c}){i,1} = relCoh';
            Cohs.std.(cond{c}){i,1} = stdCoh';
            %names.(cond{c}){i,1} = files.(cond{c}){i};
        end
        % Calculates percent change from baseline (rest) to binge
        PSDs.(cond{c}).Bands{i,3} = (PSDs.(cond{c}).Bands{i,1}-PSDs.(cond{c}).Bands{i,2})./abs(PSDs.(cond{c}).Bands{i,2});
        clear psdTrls relCoh stdCoh
    end
end
%% Define group membership
shellReduct = [-0.069767442,-0.4140625,0.06993007,-0.564343164,-0.499866986,-0.20657277,-0.4353683,-0.452188799,-0.060344828,-0.593134139,0.016771488,-0.350282486];
coreReduct = [-0.444861215,-0.15625,-0.27972028,0.105898123,-0.396212933,-0.241622575,-0.073514602,0.093333333,-0.340761374,-0.568500539,0.228710462,-0.249753208];
shellRespond = (shellReduct <= -0.26);
coreRespond = (coreReduct <= -0.26);
shellStrict = [0,1,0,1,0,0,1,1,0,0,0,0];
coreStrict = [1,0,1,0,0,0,0,0,1,0,0,1];
all = ones(length(shellReduct),1);

allGroups = {'all' 'shellRespond' 'shellStrict' 'coreRespond' 'coreStrict'};
for i = 1:length(allGroups)
    groupLog(:,i) = eval(allGroups{i})';
end

%% Extract data from each animal into vectors of all the same frequency bands from the same channel or channelpair
extract = @(C,k) cellfun(@(c)c(k),C);

for f = 1:numel(cond)
    cb.(cond{f}) = cell(numel(PSDs.(cond{f}).Bands{1}),1);
    for i = 1:numel(PSDs.(cond{f}).Bands{1})
        cb.(cond{f}){i} = extract(PSDs.(cond{f}).Bands,(i));
    end
    pb.(cond{f}) = cell(numel(Cohs.rel.(cond{f}){1}),1);
    for i = 1:numel(Cohs.rel.(cond{f}){1})
        pb.(cond{f}){i} = extract(Cohs.rel.(cond{f}),(i));
    end
end
%% Create names for each variable
% Channel-band names
bands = {'thet','alph','bet','lgam','hgam'};
channels = {'NASL','NASR','NACL','NACR'};
cbNames = cell(1,numel(cb));
for i = 1:length(channels)
    for j = 1:length(bands)
        thisName = strcat('pwr_',channels{i},'_',bands{j});
        cbNames{j+((i-1)*length(bands))} = thisName;
    end
end
% Pair-band names
chlComb = {'NASL_NASR'; 'NASL_NACL'; 'NASL_NACR'; 'NASR_NACL';'NASR_NACR';'NACL_NACR'};
pbNames = cell(1,numel(pb));
for i = 1:length(chlComb)
    for j = 1:length(bands)
        thisName = strcat('coh_',chlComb{i},'_',bands{j});
        pbNames{j+((i-1)*length(bands))} = thisName;
    end
end
%% Create data table
tic
for f = 1:numel(cond)
    % Initialize with animal name and group membership
    T.(cond{f}) = array2table(groupLog,'VariableNames',allGroups,'RowNames',{'H10','H13','H14','H15','I1','I2','I3','I4','I6','I8','I11','I12'});
    %T.(cond{f}) = table(shellRespond',shellOnly',coreRespond',coreOnly','VariableNames',{'shellRespond','shellOnly','coreRespond','coreOnly'},'RowNames',{'H10','H13','H14','H15','I1','I2','I3','I4','I6','I8','I11','I12'});
    % Add binge response columns
    T.(cond{f}) = horzcat(T.(cond{f}),table(coreReduct',shellReduct','VariableNames',{'coreReduct','shellReduct'}));
    % Fill in power data (channel-band (cb) data)
    for i = 1:numel(cbNames)
        thiscbT = table(cb.(cond{f}){i}(:,3),'VariableNames',cellstr(cbNames(i)));
        T.(cond{f}) = horzcat(T.(cond{f}),thiscbT);
    end
    % Fill in coherence data (pair-band (pb) data)
    for i = 1:numel(pbNames)
       thispbT = table(pb.(cond{f}){i},'VariableNames',cellstr(pbNames(i))); 
       T.(cond{f}) = horzcat(T.(cond{f}),thispbT);
    end
end
toc
%% Define predictor (x) and response (y) variables for regressions
allPredict = horzcat(cbNames,pbNames);
allResponse = {'shellReduct','coreReduct'};
