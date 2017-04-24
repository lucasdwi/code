function [newX] = kMean(k,x)
%% Averages each cluster of k data of x together

% INPUTS:
% k = number of datapoints to average together; format = integer
% x = row vector to be averaged

% OUTPUTS:
% newX = new, averaged, x vector

% EXAMPLE
% newx = kMean(2,[1:1:10],2)
% Will create a 5 value vector of means of every 2 data points [1.5;3.5;5.5;7.5;9.5]
%%
% Check if input is columns, if so, transpose
if iscolumn(x)
    x = x';
end

% Preallocate newX
% newX = zeros(size(x));
kCount = 0; newX = [];
for ii = 1:k:size(x,2)
    kCount = kCount+1;
    newX(kCount) = mean(x(ii:ii+k-1));
end