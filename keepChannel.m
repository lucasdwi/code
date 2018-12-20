function keepChannel(sdir,searchStr,varargin)
%% Cycles through set of files, plots data, and saves only channels given 
% by user.
% WARNING: WILL OVERWRITE FILES! Recommend making a backup in case
% something goes horrifically wrong.
%__________________________________________________________________________
% INPUTS:
% sdir = source directory to grab files from; format: string
% searchStr = string to search filenames for; format: string
% varagin = optional input, preset channels to keep, used if user knows a
%   priori that all files have the same channels to keep; format: row 
%   vector
%__________________________________________________________________________
% USE:
% rmvChannel('C:\Users\Lucas\Desktop\data\','foo')
% Will go through all files in a desktop folder 'data', with 'foo' in
% filename, plot data, and keep channels given by user.

% rmv('C:\Users\Lucas\Desktop\data\','foo',[2,4,6,8])
% Will go through all files in a desktop folder 'data', with 'foo' in
% filename, and keep channels 2,5,6, and 8.
%__________________________________________________________________________
% LLD 2017
%% Check for channels as varargin
if nargin == 3
    chans = varargin{1};
end
%% Remove channels from LFPTS.data and .label
files = [];
cd(sdir)
thisF = dir(strcat('*',searchStr,'*'));
% Just pull out 'name' field and concatenate to 'files'
files = vertcat(files,extractfield(thisF,'name')');
for fi = 1:size(files,1)
    disp(['Loading file ',num2str(fi), ' of ',num2str(size(files,1))])
    load(files{fi},'eventTs','pl2','LFPTs','adfreq')
    nC = size(LFPTs.data,1);
    figure
    for iC = 1:nC
        % Plots each channel, by default uses 2 rows and chans/2 columns
        subplot(2,nC/2,iC)
        % Plots only first 10 seconds of data for speed
        t = adfreq*10;
        plot((1:t)./adfreq,LFPTs.data(iC,1:t))
        % Constrains y-axis to ±2mV for easier differentiation of good vs.
        % bad channels
        ylim([-2 2])
        title([num2str(iC),': ',LFPTs.label{iC}])
    end
    % If channels to keep was not predefined, ask user.
    if nargin == 2
        chans = input(['Which channels do you want to keep? Input '...
        'should be in brackets, e.g. [2,4,6,8]. \nIf all channels are '...
        'to be kept, give empty set [].']);
    end
    % Only saves over file if chans to keep is given
    if ~isempty(chans)
        disp('Saving...')
        LFPTs.data = LFPTs.data(chans,:);
        LFPTs.label = LFPTs.label(1,chans);
        % Save all variables except for those used internally
        save(files{fi},'-regexp',...
        '^(?!(thisF|files|sdir|searchStr|varargin|chans|fi|iC|nC)$).')
    else
        disp('No channels selected to be saved/removed.')
    end
    close all
end