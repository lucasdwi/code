% Syntax
%   [Param, err] = PS_GetRectParam(StimN, ChN)
% Description
%   Gets parameters of the rectangular pulse for a channel ChN.
%       StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1)
%   Returns:
%       Param – array 1x5 containing parameters of the rectangular pulse: 
%           Param [0] is first phase amplitude, 
%           Param [1] is second phase amplitude, 
%           Param [2] is first phase width, 
%           Param [3] is second phase width, 
%           Param [4] is interphase delay.
%       err:
%           -1 - invalid argument(s)
%           0 - OK;
