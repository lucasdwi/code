% Syntax
%   err = PS_SetPatternType(StimN, ChN, type)
% Description
%   Configures channel ChN to use rectangular pulse parameters 
%   (pattern = 0) or preloaded arbitrary waveform pattern (pattern = 1).
%       StimN  - stimulator number to configure (starts from 1)
%       ChN  - channel number to configure (starts from 1) 
%       type - pattern type to use:
%           0 - for rectangular pulse parameter, 
%           1 - for preloaded arbitrary waveform pattern.
%   Returns:
%       -1 - invalid argument(s)
% 		0 - OK
