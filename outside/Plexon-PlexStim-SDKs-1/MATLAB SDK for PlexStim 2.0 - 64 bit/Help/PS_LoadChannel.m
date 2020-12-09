% Syntax
%   err = PS_LoadChannel(StimN, ChN)
% Description
%   Loads parameters of channel ChN to the stimulator hardware.
%   StimN  - stimulator number to configure (starts from 1)
%   ChN  - channel number which parameters will be loaded in the stimulator hardware (starts from 1)
%   Returns:
% 		-1 - invalid argument(s)
%       0 - OK
%       1 - device error
% 		3 - CRC Error
% 		6 - stimulation pattern is not ready (in case of loading an arbitrary pattern)
