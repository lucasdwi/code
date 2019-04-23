% Syntax
%   err = PS_SetRepetitions(StimN, ChN, NRep)
% Description
%   Sets number of repetitions - the number of times that the bi-phasic pulse or the arbitrary waveform (loaded from a text file) is repeated for channel ChN.
% 	Default value is 1.
%       StimN  - stimulator number to configure (starts from 1)
%       ChN  - channel number to configure (starts from 1) 
%       NRep - number of repetitions, can range from 1 to 32767; use 0 for an infinite number of repetitions.
%   Returns:
%       -1 - invalid argument(s)
%       0 - OK;
