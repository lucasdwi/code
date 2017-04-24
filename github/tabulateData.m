function [T,allPredict,allResponse,allGroups] = tabulateData(files,comp,n)
%% Collects data from group of files and puts in table (subject X variable) 
% N.B.: change 'files' structure as needed and make sure current directory
% is correct

% Inputs:
% files = array of file names to be tabulated; format = string array
%   N.B.: can be gotten from fileCycle.m
% comp = condition(s) to analyze; 1 value = from one condition; 2 values =
%   from two conditions (i.e. binge vs rest); if 2 conditions, make sure
%   there are the same number of animals (files) in each condition
% n = number of files per animal; format = integer

% Outputs:
% T = data table (subject X variable)
% allPredict = vector of variable names that can be used as predictors
% allResponse = vector of variables that are responses, predicted
% allGroups = array of strings for groups that subjects can belong to

% Example:
% [...] = tabulateData([3],{'Base'})
%   Will process all files in directory with 'Base' in name and 

% files.Base = {'H10BaseSep27','H13BaseSep27','H14BaseSep29','H15BaseSep29','I1BaseNov9','I2BaseNov9','I3BaseNov9','I4BaseSep24','I6BaseSep18','I8BaseSep23','I11BaseOct30','I12BaseNov12'};
% files.Dep = {'H10FoodDepSep25','H13FoodDepSep25','H14FoodDepSep28','H15FoodDepSep28','I1FoodDepNov13','I2FoodDep24Dec16','I3FoodDep24Nov2','I4FoodDepSep30','I6FoodDepSep30','I8FoodDepSep30','I11FoodDep24Dec16','I12FoodDepNov13'};
%% Set path
% cd(sdir);
% %% Setup files to be processed
% for f = 1:length(fType)
%     files.(fType{f}) = extractfield(dir(strcat('*',fType{f},'*')),'name')';
% end
%% Extract data from all processed files
cond = fieldnames(files);
PSDs = {};
Cohs = {};
ntrials = {};
for c = 1:length(cond)
     if length(comp) == 1
         events = {'event1'};
%         filestr = '_rest.mat';
     end
     if length(comp) == 2
         events = {'event1','event2'};
%         filestr = '_processed.mat';
     end
    for i = 1:length(files.(cond{c}))
        load(files.(cond{c})(i).name);
        for j = 1:length(events)
            PSDs.(cond{c}).Overall{i,j} = psdTrls.(events{j}).Overall;
            PSDs.(cond{c}).Bands{i,j} = psdTrls.(events{j}).Avg;
            ntrials.(cond{c}){i,j} = length(psdTrls.(events{j}).Pow);
            Cohs.rel.(cond{c}){i,1} = coh.relCoh';
            if length(comp) == 1
                PSDs.(cond{c}).Bands{i,length(events)+1} = relPower;
            end
            if length(comp) == 2
                Cohs.std.(cond{c}){i,1} = stdCoh';
            end
        end
        % Calculates percent change from baseline (rest) to binge
        if length(comp) == 2
            PSDs.(cond{c}).Bands{i,3} = (PSDs.(cond{c}).Bands{i,1}-PSDs.(cond{c}).Bands{i,2})./abs(PSDs.(cond{c}).Bands{i,2});
        end
        clear psdTrls relPower relCoh stdCoh;
    end
end
%% Define group membership
% Order: H10,H13,H14,H15,I11,I12,I1,I2,I3,I4,I6,I8
shellPercents = [-0.069767442,-0.4140625,0.06993007,-0.564343164,0.016771488,-0.350282486,-0.499866986,-0.20657277,-0.4353683,-0.452188799,-0.060344828,-0.593134139];
corePercents= [-0.444861215,-0.15625,-0.27972028,0.105898123,0.228710462,-0.249753208,-0.396212933,-0.241622575,-0.073514602,0.093333333,-0.340761374,-0.568500539];
% Replicate Percents for each n files per animal
shellReduct = repmat(shellPercents,n,1);
shellReduct = shellReduct(:)';
coreReduct = repmat(corePercents,n,1);
coreReduct = coreReduct(:)';
% Setup responder group based on a reduction of more than 26%
shellRespond = (shellReduct <= -0.26);
coreRespond = (coreReduct <= -0.26);
% Replicate Stricts for each n files per animal
shellStrict = repmat([0,1,0,1,0,0,0,0,1,1,0,0],n,1);
shellStrict = shellStrict(:)';
coreStrict = repmat([1,0,1,0,0,1,0,0,0,0,1,0],n,1);
coreStrict = coreStrict(:)';

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
for f = 1:numel(cond)
    % Extract file names from files structure
    fNames = extractfield(files.(cond{f}),'name')';
    % Initialize with animal name and group membership
    T.(cond{f}) = array2table(groupLog,'VariableNames',allGroups,'RowNames',fNames);
    %T.(cond{f}) = table(shellRespond',shellOnly',coreRespond',coreOnly','VariableNames',{'shellRespond','shellOnly','coreRespond','coreOnly'},'RowNames',{'H10','H13','H14','H15','I1','I2','I3','I4','I6','I8','I11','I12'});
    % Add binge response columns
    T.(cond{f}) = horzcat(T.(cond{f}),table(coreReduct',shellReduct','VariableNames',{'coreReduct','shellReduct'}));
    % Fill in power data (channel-band (cb) data)
    for i = 1:numel(cbNames)
        if length(comp) == 1
            thiscbT = table(cb.(cond{f}){i}(:,2),'VariableNames',cellstr(cbNames(i)));
        end
        if length(comp) == 2
            thiscbT = table(cb.(cond{f}){i}(:,3),'VariableNames',cellstr(cbNames(i)));
        end
        T.(cond{f}) = horzcat(T.(cond{f}),thiscbT);
    end
    % Fill in coherence data (pair-band (pb) data)
    for i = 1:numel(pbNames)
       thispbT = table(pb.(cond{f}){i},'VariableNames',cellstr(pbNames(i))); 
       T.(cond{f}) = horzcat(T.(cond{f}),thispbT);
    end
end
%% Define predictor (x) and response (y) variables for regressions
allPredict = horzcat(cbNames,pbNames);
allResponse = {'shellReduct','coreReduct'};

