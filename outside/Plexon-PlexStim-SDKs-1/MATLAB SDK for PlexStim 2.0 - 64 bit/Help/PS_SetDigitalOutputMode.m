% Syntax
%   err = PS_SetDigitalOutputMode(StimN, Mode)
% Description
%   Sets Digital Output mode to low or high during inter-pulse interval.
%   Each stimulator channel has a dedicated digital output that indicates when stimulation is occurring on that channel. 
%   The digital output is always high during the pulse or arbitrary waveform output, but the user can control the state 
%   of the digital output during the time in between pulses or arbitrary waveforms by using this function.
% 	Digital output can be set to low or high during inter-pulse interval.
%   Default value is 1 (low).
%       StimN  - stimulator number to configure (starts from 1)
%       Mode – can have the following values:
%           0 – high;
%           1 – low.
%       Returns:
%           -1 - invalid argument(s)
%           0 - OK;
   