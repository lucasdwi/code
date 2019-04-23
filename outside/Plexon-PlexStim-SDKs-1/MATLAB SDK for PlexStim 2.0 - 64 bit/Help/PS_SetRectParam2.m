% Syntax
%   err = PS_SetRectParam2(StimN, ChN, Param)
% Description
%   Sets parameters of the rectangular pulse for a channel. Has the same effect as PS_SetRectParam, but uses a structure as the third input argument instead of an array.
%       StimN  - stimulator number to configure (starts from 1)
%       ChN  - channel number to configure (starts from 1)
%   Param – structure containing parameters of the rectangular pulse; the structure has 5 fields: 
%       Pattern.A1 is first phase amplitude, 
%       Pattern.A2 is second phase amplitude, 
%       Pattern.W1 is first phase width,
%       Pattern.W2 is second phase width,
%       Pattern.Delay is interphase delay.
%   The structure should be defined, for example, as  
%       Pattern.A1 = 100;
%       Pattern.A2 = -100;
%       Pattern.W1 = 25;
%       Pattern.W2 = 25;
%       Pattern.Delay = 25;
% 
%   before calling PS_SetRectParam2 function.
% 
%   The default values are:
%       - first phase amplitude = 100 mA;
%       - second phase amplitude = -100 mA;
%       - first phase width = 50 us;
%       - second phase width = 50 us;
%       - interphase delay = 25 us;
%   Returns:
%    	-1 - invalid argument(s)
% 		0 - OK;
