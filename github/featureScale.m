function [scaledX] = featureScale(x,minX,maxX)
%% Scales data x from 0 to 1
% INPUTS
% x = data to be scaled; format = matrix with data in rows belonging to
%   same set
% varargin{1} = min
% varargin{2} = max
%%
% Use min and max from x
if isempty(minX) && isempty(maxX)
    minX = min(X); maxX = max(x);
else
% Use given min and max
scaledX = (x-minX)/(maxX-minX);