function [x,tvec] = sigGen(fs,t,f,a,phi)
%% Generates either one signal in the case of unit inputs, or combines many oscillations if given vectors
% INPUTS
% fs = sampling rate
% t = length of time to use
% f = frequency of signal; if vector applies each value to its own
%       signal
% a = amplitude of signal; if vector applies each value to its own
%       signal
% phi = phase offset of signal in degrees; if vector applies each value to
%       its own signal
% OUTPUTS
% x = signal
% t = time vector
%% Check amplitudes and phase offsets
% Check if a is empty, if so set it to ones of equal size as f
if isempty(a)
    a =  ones(size(f,2));
    % If a is not empty but smaller than f repeat to size f
elseif length(a) < size(f,2)
    a = repmat(a,1,size(f,2));
end
% If phi is empty, set it to zeros of equal size as f
if isempty(phi)
    phi = zeros(size(f,2));
elseif size(phi,2) < size(f,2)
    phi = repmat(phi,1,size(f,2));
end
% Convert phi degrees to radians
phi = (phi*pi)/180;
%%
% Create time vector from fs and t
tvec = 0:1/fs:t;
xMat = zeros(length(f),t*fs+1);
for fi = 1:length(f)
    xMat(fi,:) = sin(2*pi*f(fi)*tvec+phi(fi))*a(fi);
end
% Combine signals
x = sum(xMat,1);