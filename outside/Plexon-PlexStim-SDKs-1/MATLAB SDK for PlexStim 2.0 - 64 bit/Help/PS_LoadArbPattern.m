% Syntax
%   err = PS_LoadArbPattern(StimN, ChN, s)
% Description
%   Load an arbitrary waveform pattern from a (.pat) file into a selected channel.
%   StimN  - stimulator number to configure (starts from 1)
%   ChN  - channel number to configure (starts from 1)
%   s – string, contains full path (not more than 512 characters) of the file with arbitrary waveform pattern
%    Returns:
%       -1 - invalid argument(s)
%       0 - pattern loaded from file successfully
%       7 - length of the file name exceeds 512
%       8 - file doesn't exist
%       9 - file is opened by another process
%       10 - number of points in the pattern exceeds maximum allowed size of 1000
%       11 - file contains invalid value(s)
%       12 – file contains too few lines (min is 3)
%       13 – file contains a mismatched amplitude duration pair











