%% Addpaths
% Code
addpath(genpath('C:\Users\GreenLab\Documents\GitHub'))
% Remove fieldtrip
% rmpath(genpath('C:\Users\Pythia\Documents\GreenLab\code\outside\fieldtrip-master\'))
% Data
% addpath(genpath('C:\Users\Pythia\Documents\GreenLab\data\'))
cd('C:\Users\GreenLab\desktop')
%% Set up plotting defaults
% Reset defaults (perhaps redundant, but better safe than sorry)
reset(0);
% Set figure fonts to Arial rather than Helvetica
set(0,'defaultTextFontname','Arial')
set(0,'defaultAxesFontName','Arial')
set(0,'defaultLegendFontName','Arial')
% Set font sizes - N.B. for multipliers of the default (12 pt
% font = 1) use quotient of desired font size by 12, e.g. font size of 18 =
% 18/12
set(0,'defaultAxesFontSize',12)
set(0,'defaultAxesTitleFontSizeMultiplier',18/12)
set(0,'defaultAxesLabelFontSizeMultiplier',14/12)
set(0,'defaultLegendFontSize',12)
set(0,'defaultTextFontSize',12)
% Set title font weights to regular rather than bold
set(0,'defaultAxesTitleFontWeight','normal')
% Set axes line and font color to black
set(0,'defaultAxesXColor','k')
set(0,'defaultAxesYColor','k')
set(0,'defaultAxesZColor','k')
% Turn off legend box
set(0,'defaultLegendBox','off')
% Set line width to 2
set(0,'defaultLineLineWidth',2)