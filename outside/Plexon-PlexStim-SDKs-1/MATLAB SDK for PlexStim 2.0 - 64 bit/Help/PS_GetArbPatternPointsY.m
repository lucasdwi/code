% Syntax
%   [CoordsY, err] = PS_GetArbPatternPointsY(StimN, ChN)
% Description
%   Gets Y coordinates of a graphical representation of the arbitrary waveform contained in the loaded pattern file.
%       StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1)
% 	Returns:
% 		CoordsY - array of Y coordinates of the points; 
% 				it accomodates NPoints of integer values;
% 				CoordsY contains coordinates as y1 y2 y3 ....
% 		err:
% 			-1 - invalid argument(s)
% 			0 - OK
%           5 - NPoints is not equal number of point for the pattern in this channel
