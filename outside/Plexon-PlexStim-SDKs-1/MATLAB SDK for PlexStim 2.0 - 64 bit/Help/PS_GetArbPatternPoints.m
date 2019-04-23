% Syntax
%   [Coords, err] = PS_GetArbPatternPoints(StimN, ChN)
% Description
%   Gets X and Y coordinates of a graphical representation of the arbitrary waveform pattern loaded into the selected stimulator and channel.
%   These coordinates can be used to draw a graph of the pattern
%       StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1)
% 	Returns:
% 		Coords - array of integer coordinates of the points; 
% 			it accomodates NPoints*2 integer values;
% 			“Coords” contains coordinates as x1 y1 x2 y2 ....
% 		err:
% 			-1 - invalid argument(s)
% 			0 - OK
%           5 - NPoints is not equal number of point for the pattern in this channel
