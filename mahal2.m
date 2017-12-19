function [dM] = mahal2(x,y)
%%
% Get sample sizes
nx = size(x,1);
ny = size(y,1);
% Get pooled sample size
N = nx + ny;
% Get covariance of both samples
Cx = cov(x);
Cy = cov(y);
% Calculate ubiased pooled covariance matrix
S = (nx*Cx + ny*Cy)/(N-2);
% Get mean of each dimension
mx = mean(x,1);
my = mean(y,1);
% Calculate Mahalanobis distance
dM = sqrt((mx-my)*S^-1*(mx-my)');