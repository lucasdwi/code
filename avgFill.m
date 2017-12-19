function [xFill,yFill] = avgFill(x,y,dim,varargin)
%% Obtains average of x and y matrices and provides new x and y vectors for
% use with the plotting function 'fill'. New xFill and yFill are calculated
% as the average plus/minus one standard deviation.
%__________________________________________________________________________
% INPUTS:
% x = matrix of x values
% y = matrix of y values
% dim = dimension to average over
% numS = optional; number of standard deviations; default = 1
%__________________________________________________________________________
% OUTPUTS:
% xFill = x vector
% yFill = y vector 
%__________________________________________________________________________
% LLD 2017
%% Check varargin
if nargin > 3
    numS = varargin{1};
else 
    numS = 1;
end
% Get mean of both x and y, along dimension given
mX = mean(x,dim);
mY = mean(y,dim);
% Get standard deviation of both x and y, along dimension given; multiply
% by number of standard deviations
sX = numS*std(x,[],dim);
sY = numS*std(y,[],dim);
% Use trig to find point that is sX and sY away from each point at 135
% degrees - i.e. that forms the upper left boundary (1 standard deviation
% below x and above y)
theta = (135*pi)/180;
x1 = mX+cos(theta).*sX;
y1 = mY+sin(theta).*sY;
% Do the same for the lower boundary using 315 degrees (1 standard devation
% above x and below y)
theta = (315*pi)/180;
x2 = mX+cos(theta).*sX;
y2 = mY+sin(theta).*sY;
% Concatenate pluses with flipped minuses
xFill = [x1;flipud(x2)];
yFill = [y1;flipud(y2)];