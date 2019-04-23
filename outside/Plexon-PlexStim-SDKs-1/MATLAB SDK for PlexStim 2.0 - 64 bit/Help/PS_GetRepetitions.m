% Syntax
%   [NRep, err] = PS_GetRepetitions(StimN, ChN)
% Description
%   Gets number of repetitions - the number of times that the bi-phasic pulse or the arbitrary waveform (loaded from a text file) is repeated for channel ChN.
% 	StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1) 
% 	Returns:
%       NRep - number of repetitions, can range from 1 to 32767; 0 for an infinite number of repetitions.
%   err:
%    	-1 - invalid argument(s)
% 		0 - OK;
