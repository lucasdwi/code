function plotMIs(MI,fs,fSpaceR,interval,rotX)

%input:

%MI is a comodulogram

%fs is the sampling rate

%fSpaceR is the frequency space (rad/sample) returned by Lilly's code

%Interval is used to separate adjacent frequencies.  Larger values generate
%courser frequency axes.

%rotX:  Set to 1 to rotate xlabels 45 deg.  Can improve readibility

%Example: Plot the comodulogram from channel 1 to channel 3 in a multichannel
%data set.  Assume data sampled at 500 Hz, we skip every other frequency in
%building the frequency axes, and we do not rotate the labels.

%plotMIs(squeeze(MI(:,:,1,3)),fs,fSpace,2,0);

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

%A. Nakhnikian, 2014

fSpace = fSpaceR*(fs/(2*pi));
MI = rot90(MI,2);
figure, hold on, contourf(MI','Linecolor','none')
plotInds = 1:interval:length(fSpace);
set(gca,'ytick',plotInds), set(gca,'xtick',plotInds)
labels = cell(size(plotInds));
for ind = 1:length(plotInds)
    labels{ind} = num2str(fSpace(plotInds(ind)),'%3.1f');
end
set(gca,'xticklabel',fliplr(labels))
set(gca,'yticklabel',fliplr(labels))
set(gca,'ydir','normal')
set(gca,'xlim',[1 length(MI)])
set(gca,'ylim',[1 length(MI)])
colorbar
set(gca,'FontWeight','bold')
xlabel('Phase Frequency (Hz)','FontWeight','bold')
ylabel('Amplitude Frequency (Hz)','FontWeight','bold')
if rotX==1
rotateXLabels(gca,45);
end
 %set(gca,'xticklabel',round(fliplr(fSpace([1:interval:length(fSpace)]))'))
% set(gca,'yticklabel',round((fSpace(1:interval:length(fSpace)))))