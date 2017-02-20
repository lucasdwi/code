function [wAvg] = wMean(x,w,dim)
%%
% INPUTS
% x = data matrix with format at least 2 dimensions
% w = vector of weights corresponding to dimension of x; dim
% dim = dimension to average over, 'collapse'
%%
% Replace any inf with NaN
x(isinf(x)) = NaN;
% Repeat w vector into 2d matrix
wMat = repmat(w',1,size(x,2));
% Stack 2dW and permute to create 3D weight array equivalent to data
wArr = permute(repmat(wMat,1,1,size(x,1)),[3,2,1]);
% Multiply data by weights
wDat = x.*wArr;
% Sum across dim
datSum = nansum(wDat,dim);
% Get sum of weights accounting for any NaNs
wArr(isnan(x)) = 0;
wSum = sum(wArr,3);
% Get weighted average
wAvg = datSum./wSum;
