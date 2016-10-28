function [newX] = kMean(k,x)
%% Averages i:k of x together

% INPUTS:
% k = number of datapoints to average together; format = integer
% x = vector of data to be averaged

% OUTPUTS:
% newX = new, averaged, x vector

% EXAMPLE
% newx = kMean(2,[1:1:10])
% Will create a 5 value vector of means of every 2 data points []
%%
for ii = 1
    for j = 1:2:24
        new(j,ii) = mean([test(j,ii),test(j+1,ii)]);
    end
end
new = new(1:2:end,:);