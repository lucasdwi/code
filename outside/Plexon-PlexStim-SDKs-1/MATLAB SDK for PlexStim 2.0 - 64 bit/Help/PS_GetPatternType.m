% Syntax
%   [type, err] = PS_GetPatternType(StimN, ChN)
% Description
%   Checks if channel set up to use rectangular pulse parameters or preloaded arbitrary waveform pattern
%   StimN  - stimulator number to query (starts from 1)
%   ChN  - channel number to configure (starts from 1) 
%   Returns:
%       -	type - which pattern type is in use: 0 if rectangular pulse, 1 if preloaded arbitrary waveform.
%       - err:
%           - 1 - invalid argument(s)
%           0 - OK.
