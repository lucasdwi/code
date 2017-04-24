function [models,Tstats,sigModels] = lineReg(T,allPredict,allResponse,allGroups,a)
%% Run linear regression on all variables with groups paired with equivalent responses (e.g. shellStrict with shellReduct but not coreReduct)
% For code for all combinations check LineRegNotes
cond = fieldnames(T);
for f = 1:numel(cond)
    mI = 0;
    for i = 1:length(allGroups)
        for k = 1:length(allPredict)
            if i == 1
                for j = 1:length(allResponse)
                    thismd = fitlm(T.(cond{f})(T.(cond{f}).(allGroups{i}) == 1,:),strcat(allResponse{j},'~',allPredict{k}));
                    mI = mI + 1;
                    thisind = mI;
                    models.(cond{f}){1,thisind} = thismd;
                    models.(cond{f}){2,thisind} = allGroups{i};
                    models.(cond{f}){3,thisind} = allPredict{k};
                    models.(cond{f}){4,thisind} = allResponse{j};
                end
            else
                if i == 2 || i == 3
                    thismd = fitlm(T.(cond{f})(T.(cond{f}).(allGroups{i})== 1,:),strcat(allResponse{1},'~',allPredict{k}));
                    mI = mI + 1;
                    thisind = mI;
                    models.(cond{f}){1,thisind} = thismd;
                    models.(cond{f}){2,thisind} = allGroups{i};
                    models.(cond{f}){3,thisind} = allPredict{k};
                    models.(cond{f}){4,thisind} = allResponse{1};
                else
                    if i == 4 || i == 5
                        thismd = fitlm(T.(cond{f})(T.(cond{f}).(allGroups{i})== 1,:),strcat(allResponse{2},'~',allPredict{k}));
                        mI = mI + 1;
                        thisind = mI;
                        models.(cond{f}){1,thisind} = thismd;
                        models.(cond{f}){2,thisind} = allGroups{i};
                        models.(cond{f}){3,thisind} = allPredict{k};
                        models.(cond{f}){4,thisind} = allResponse{2};
                    end
                end
            end
        end
    end
end
%% Create stats table
% Setup empty table
for f = 1:numel(cond)
    Tstats.(cond{f}) = cell2table({});
    for i = 1:length(models.(cond{f}))
        predictStr = [];
        for j = 1:length(models.(cond{f}){1,i}.PredictorNames)
            if j == 1
            thisPredict = strcat(models.(cond{f}){1,i}.PredictorNames(j));    
            else
            thisPredict = strcat('+',models.(cond{f}){1,i}.PredictorNames(j));
            end
            predictStr = strcat(predictStr,thisPredict);
        end
        thisName = strcat(models.(cond{f}){2,i},'_',models.(cond{f}){1,i}.ResponseName,'~',predictStr);
        thisTstat = table(models.(cond{f}){1,i}.Rsquared.Ordinary',models.(cond{f}){1,i}.Rsquared.Adjusted',coefTest(models.(cond{f}){1,i})',i,{models.(cond{f}){2,i}},'RowNames',thisName,'VariableNames',{'R2_Ord','R2_Adj','pValue','ModelNumber','Group'});
        Tstats.(cond{f})(thisName,:) = thisTstat;
    end
end
%% Find models with significant p values
for f = 1:numel(cond)
    sigModels.(cond{f}) = Tstats.(cond{f})(Tstats.(cond{f}).pValue <= a,:);
    sigModels.(cond{f}) = sortrows(sigModels.(cond{f}),'pValue','ascend');
end
