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
% Get mX +/- sX
xP = mX+sX;
xM = mX-sX;
% Get mY +/- sY
yP = mY+sY;
yM = mY-sY;
% Concatenate pluses with flipped minuses
xFill = [xP;flipud(xM)];
yFill = [yP;flipud(yM)];