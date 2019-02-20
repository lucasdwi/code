% Syntax
%   err = PS_SetRectParam(StimN, ChN, Param)
% Description
%   Sets parameters of the rectangular pulse for a channel
%       StimN  - stimulator number to configure (starts from 1)
%       ChN  - channel number to configure (starts from 1)
%       Param – array 1x5 containing parameters of the rectangular pulse; it should be defined, for example, as  
%       pattern = [100, -100, 25, 25, 25] before calling PS_SetRectParam function.
% 
%       Param [0] is first phase amplitude, 
%       Param [1] is second phase amplitude, 
%       Param [2] is first phase width, 
%       Param [3] is second phase width, 
%       Param [4] is interphase delay.
%   The default values are:
%       - first phase amplitude = 100 mA;
%       - second phase amplitude = -100 mA;
%       - first phase width = 50 us;
%       - second phase width = 50 us;
%       - interphase delay = 25 us;
%   Returns:
%    	-1 - invalid argument(s)
% 		0 - OK;
