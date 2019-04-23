% Syntax
%   err = PS_SetPeriod(StimN, ChN, Period)
% Description
%   Set period for a channel in milliseconds
% 	Default value is 5 ms.
%       StimN  - stimulator number to configure (starts from 1)
%       ChN  - channel number to configure (starts from 1) 
%       Period - period value in milliseconds, valid values are from 0.020 ms <= Period <= 125,000 ms
%   Returns:
%       -1 - invalid argument(s)
%       0 - OK;
