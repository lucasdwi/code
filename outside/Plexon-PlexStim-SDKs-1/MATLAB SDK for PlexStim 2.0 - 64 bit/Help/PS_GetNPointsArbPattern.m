% Syntax
%   [NPoints, err] = PS_GetNPointsArbPattern(StimN, ChN)
% Description
%   Get the number of points in a graphical representation of the arbitrary waveform.
%       StimN  - stimulator number to query (starts from 1)
%       ChN  - channel number to query (starts from 1)
%   Returns:
%       NPoints - the number of points in a graphical representation of the arbitrary waveform pattern loaded into the selected stimulator and channel
%       err:
%           -1 - invalid argument(s)
%           0 - OK.
