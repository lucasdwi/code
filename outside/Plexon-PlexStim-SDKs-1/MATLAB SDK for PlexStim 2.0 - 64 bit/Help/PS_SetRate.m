% Syntax
%   err = PS_SetRate(StimN, ChN, Rate)
% Description
%   Set repetition rate for a channel in Hertz.
% 	Default value is 200 Hz.
%       StimN  - stimulator number to configure (starts from 1)
%       ChN  - channel number to configure (starts from 1) 
%       Rate - rate in Hertz; valid values are from 0.008 Hz <= Rate <= 50000 Hz
%   Returns:
%       -1 - invalid argument(s)
%       0 - OK;
