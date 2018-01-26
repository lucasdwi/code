function [sdir,file,filter,dsf,thresh,onset,offset,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent] = scbParamsSingle
%% Populates workspace with input variables for spectcompbase.m
%__________________________________________________________________________
% INPUTS:
% Altered in function itself.
%__________________________________________________________________________
% OUTPUTS:
% sdir = source directory of data; format: string file = file name w/o
% extension; format: string filter = wether or not to use a filter; format:
% 'y' or 'n' dsf = factor with which to downfactor; format: whole integer
% thresh = threshold for detecting noise artifacts; format: mV (or other
%   y-axis scale for time series)
% onset = number of samples to NaN before noise event; format: number
%  in discrete samples, e.g. samples = sec * adfreq
% offset = number of samples to NaN after noise event; format: number
%   in discrete samples, e.g. samples = sec * adfreq
% minInt = minimum interval of clean data to use, N.B.: consider the lowest
%   frequency of interest and coherence method; format: seconds
% foi = frequencies of interest; format = [lower step upper] in Hz bands =
% structure with bands of interest starting with lowest; format:
%   {'band1',[lower upper];'band2',[lower upper];...}
% cycles = number of cycles to use in creating frequency dependent windows
%   in ft_freqanalysis; format: number
% ftimwin = size of window to use in ft_freqanalysis; format: number in
%   seconds; N.B.: during PowerCorr.m ftimwin will be used to compute the
%   number of cycles at the lowest frequency band of interest and will warn
%   if less than 3
% overlap = amount of overlap to use with sliding windows; format: percent
%   in decimal form (1-percent; e.g. 90% overlap = 0.1)
% eoi = events of interest, if one event then normalizes within that event,
%   if two then compares events; format: structure {'tag1',[0 3];'tag2',[0
%   3]} first column corresponds to event marker labels, second column
%   corresponds to timing around given event. If one number is given then
%   it will be treated as the minimum trial length and as many trials of
%   that length as possible will be found. If two numbers are given then
%   the 'window' around the given maker will be found. N.B.: if all the
%   data use the tag 'full', otherwise use tags corresponding to event
%   markers: 'app','binge','rest','sleepOut'/'sleepIn'
% saveParent = parent directory path to save plots and files; format:
%   string N.B.: if directory/file exists will warn about overwritting
%__________________________________________________________________________
%% LLD 2017
%% See above documentation for details on variables.

sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\angela\mat\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\raw\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\mike\toProcess\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\baseData\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\';
%sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\editedWithSleep\';
% file = 'N07_4_21_16_Base_pl2_scored_plx';
% file = 'M376_Dec_6_Base_ScoredComplete';
file = 'MC8_2017-06-01_10s10r_Anesth.mat';
% file = '10HzR_N12_PostStim_9_20_16';
filter = 'y';
dsf = 5;
thresh = 2; 
onset = 0.0125;
offset = 5;
% minInt = 5;
foi = [1 2 100];
% foi = [70 2 90];
% bands = {'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
% bands = {'hgam',[70 90]};
cycles = 3;
ftimwin = [];
overlap = 0.5;
cohMethod = 'mat';
eoi = {'Base',[0 5];'Int1',[0 5];'Int2',[0 5];'Int3',[0 5];'Int4',[0 5];'Int5',[0 5];'Post',[0 5]};
% eoi = {'binge',[0 5];'notbinge',[0 5]};
% eoi = {'binge (s',[-5 0]};%'binge (s',[-6 -1];'binge (s',[-7 -2];'binge (s',[-8 -3];'binge (s',[-9 -4];'binge (s',[-10 -5]};
saveParent = 'C:\Users\Lucas\Desktop\GreenLab\data\paper2\test\';
