function [newX] = normOversample(x,n) 
%% Fits a normal distribution to each column of x and then randomally 
% samples that distribution n times
%__________________________________________________________________________
% INPUTS:
% x = data matrix; format = observation X variable
% n = number of samples to generate; format = integer
%__________________________________________________________________________
% OUTPUTS:
% newX = matrix of oversampled data; format = n X original number of
%   variables
%__________________________________________________________________________
% LLD 2017
%% Rand dist and imputer (manual bootstrapping)
newX = zeros(n,size(x,2));
for ii = 1:size(x,2)
   m = mean(x(:,ii));
   s = std(x(:,ii));
   pd = makedist('Normal','mu',m,'sigma',s);
   newX(:,ii) = random(pd,1,n)';
end
