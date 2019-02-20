% Syntax
%   [Value, err] = PS_IsWaveformBalanced(StimN, ChN)
% Description
%   Checks if the stimulation waveform is balanced. 
%   Both rectangular pulses and arbitrary patterns should be charge balanced. It means that the same amount of current is deposited and withdrawn from the electrode. This function analyzes the pulse or pattern to check the net charge is not zero.Gets duration of the whole stimulation pattern.
%       StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1) 
%   Returns:
%       Value: 1 if the waveform is balanced (net charge is zero), 0 otherwise
%   err:
%       -1 - invalid argument(s)
% 		0 - OK;
