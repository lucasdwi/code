function [T,allPredict,allResponse,models,Tstats,sigModels] = regFlow(a)
%% First tabulate data as defined within tabulateData.m
[T,allPredict,allResponse,allGroups] = tabulateData;
%% Then run linear regression on table T with allVars being predictors
[models,Tstats,sigModels] = lineReg(T,allPredict,allResponse,allGroups,a);