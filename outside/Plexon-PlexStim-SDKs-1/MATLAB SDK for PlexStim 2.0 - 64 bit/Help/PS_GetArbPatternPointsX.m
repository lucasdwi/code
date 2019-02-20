% Syntax
%   [CoordsX, err] = PS_GetArbPatternPointsX(StimN, ChN)
% Description
%   Gets X coordinates of a graphical representation of the arbitrary waveform contained in the loaded pattern file.
%       StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1)
% 	Returns:
% 		CoordsX - array of X coordinates of the points; 
% 			it accomodates NPoints of integer values;
% 			CoordsX contains coordinates as x1 x2 x3 ....
% 		err:
% 			-1 - invalid argument(s)
% 			0 - OK
%           5 - NPoints is not equal number of point for the pattern in this channel
