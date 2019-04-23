% Syntax
%   [scaling, err] = PS_GetVmonScaling(StimN)
% Description
%   Gets the Voltage monitor scaling for stimulation mode.
% 	StimN  - stimulator number to monitor (starts from 1)
%   Returns:
%       -Scaling - scaling parameter (0.25 V/V, 2.5 V/V, 25 V/V, or 250 V/V).
%       -1 - invalid argument(s)
% 		0 - OK
