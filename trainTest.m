function [trainX,trainY,testX,testY,trainInd,testInd] = trainTest(x,y,perc)
%% Splits dataset into training and testing sets
% INPUTS
% x = x data to be split; format: array, samples X features
% y = y data to be split; format: array
% perc = percent of data to put into test set; format: decimal (0.2 = 20%)
% OUTPUTS
% trainX = x data for training
% trainY = y data for training
% testX = x data for testing
% testY = y data for testing
%__________________________________________________________________________
% LDD 2019-12-16
%%
n = size(x,1);
trainN = floor(n*(1-perc));
testN = n-trainN;
inds = randperm(n,n);
trainInd = inds(1:trainN);
testInd = inds(trainN+1:end);
trainX = x(trainInd,:);
trainY = y(trainInd);
testX = x(testInd,:);
testY = y(testInd);