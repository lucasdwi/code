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
% Calculate size of input data
n = size(x,1);
% Calculate number of samples for training
trainN = floor(n*(1-perc));
% Calculate number of samples for testing
testN = n-trainN; % not used at this point
% Create randomized index list 
inds = randperm(n,n);
% Split randomized index list into train and test indices 
trainInd = inds(1:trainN);
testInd = inds(trainN+1:end);
% Use train and test indices to split data into train and test sets
trainX = x(trainInd,:);
trainY = y(trainInd);
testX = x(testInd,:);
testY = y(testInd);