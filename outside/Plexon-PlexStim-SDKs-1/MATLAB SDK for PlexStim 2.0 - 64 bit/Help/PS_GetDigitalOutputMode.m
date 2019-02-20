% Syntax
% 	[Mode, err] = PS_GetDigitalOutputMode(StimN)
% Description
% 	Checks if the Digital Output is low or high during the inter-pulse
% 	interval.  Each stimulator channel has a dedicated digital output that
% 	indicates when stimulation is occurring on that channel. The digital
% 	output is always high during the pulse or arbitrary waveform output,
% 	but the user can control the state of the digital output during the
% 	time in between pulses or arbitrary waveforms.
% 	Default value is 1 (low).
% 
% 	StimN  - stimulator number to query (starts from 1)
% 
% 	Returns:
% 		Mode - can have the following values:
% 			0 – high;
% 			1 – low.
% 		err:
% 			-1 - invalid argument(s)
% 			0 - OK;


















