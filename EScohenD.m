function [d] = EScohenD(data1,data2)
%% Caclulate Cohen's D effect size using pooled standard deviation (Hedges' g)
% INPUTS:
% data1 = dataset 1; format = row or column vector
% data2 = dataset 2; format = row or column vector

% OUTPUT:
% d = Cohen's D
%%
% Get basic stats from datasets
mu1 = mean(data1); mu2 = mean(data2);
s1 = std(data1); s2 = std(data2);
n1 = length(data1); n2 = length(data2);
% Calculate pooled standard deviation (Hedges' g)
s = sqrt((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2);
% Calculate Cohen's D
d = abs((mu1 - mu2))/s;