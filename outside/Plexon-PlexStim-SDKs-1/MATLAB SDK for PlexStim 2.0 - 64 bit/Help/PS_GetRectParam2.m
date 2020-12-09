% Syntax
%   [Param, err] = PS_GetRectParam2(StimN, ChN)
% Description
%   Gets parameters of the rectangular pulse for a channel ChN.
%   Has the same effect as the function PS_SetRectParam2, but returns parameters of the pul;se as a structure having 5 fields.
%       StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1)
%   Returns:
%       Param: structure having 5 fields:
%       Pattern.A1 is first phase amplitude, 
%       Pattern.A2 is second phase amplitude, 
%       Pattern.W1 is first phase width,
%       Pattern.W2 is second phase width,
%       Pattern.Delay is interphase delay.
% 	err:
%    	-1 - invalid argument(s)
% 		0 - OK;
