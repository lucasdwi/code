function [amplitudeSeries,phaseSeries] = waveletBandPass(x,y,psi)

%FUNCTION: Returns instantaneous amplitude (for x) and phase (for y) of
%each time series.  Generates analytic time series of x and y at each
%frequency under analysis using the wavelet(s) contained in the array psi.
%Psi is generated using morsespace and morsewave from the Jlab analysis
%package.

%License:

%  This software is distributed under the "Creative Commons Attribution
%  Noncommercial-Share Alike License"
%
%     Version 3.0, available at
%
%         http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
%     You are free:
%
%         To Share -- To copy, distribute and transmit the work
%
%         To Remix -- To adapt the work
%
%     Under the following conditions:
%
%         Attribution -- You must attribute the work in the manner specified
%            by the author or licensor (but not in any way that suggests that
%            they endorse you or your use of the work).
%
%         Noncommercial -- You may not use this work for commercial purposes.
%
%         Share Alike -- If you alter, transform, or build upon this work,
%            you may distribute the resulting work only under the same or
%            similar license to this one.
%
%     See the above link for the full text of the license.
%     _______________________________________________________________________
%
%     Disclaimer
%
%     This software is provided 'as-is', without any express or implied
%     warranty. In no event will the author be held liable for any damages
%     arising from the use of this software.

%A. Nakhnikian Feb 2014

if size(x,2)>size(x,1)
    x = x';
end
if size(y,2)>size(y,1)
    y = y';
end
X = wavetrans(x,psi);
Y = wavetrans(y,psi);
% if numel(size(X)) == 4||numel(size(x))<3&&numel(size(psi))>2
%     %make each eigenspectrum a "trial"
%     X = reshape(X,size(X,1),size(X,2),size(X,3)*size(X,4));
%     Y = reshape(Y,size(Y,1), size(Y,2), size(Y,3)*size(Y,4));
%     %Y = repmat(Y(:,:,:,1),1,1,size(psi,3));
% end
amplitudeSeries = abs(X);
phaseSeries = angle(Y);
end