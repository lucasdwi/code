function [pow,value] = nearestPow2(x)
%% Finds the nearest power of two to x; if want the next largest use nextpow2
% INPUTS:
% x = number to find the nearest power of two of

% OUTPUTS:
% pow = nearest power of two
% value = value; 2^power
%%
pUp = nextpow2(x);
pDown = pUp-1;
p = [pUp,pDown];
[~,ind] = min([abs(x-2^pUp),abs(x-2^pDown)]);
pow = p(ind);
value = 2^pow;