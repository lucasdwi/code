% Syntax
%   [Value, err] = PS_GetStimPatternDuration(StimN, ChN)
% Description
%   Gets duration of the whole stimulation pattern.
%   StimN  - stimulator number to query (starts from 1)
%   ChN  - channel number to query (starts from 1) 
%   Returns:
%       Value - duration of the stimulation pattern in microseconds;
%   err:
%       -1 - invalid argument(s)
% 		0 - OK;
