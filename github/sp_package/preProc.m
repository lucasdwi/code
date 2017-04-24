function [] = spectcompbase(sdir,file,filter,dsf,thresh,onset,offset,minInt,foi,bands,eoi)
% ('C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\','N7_PreStim_9_8_16','y',5,2.5,5,17000,3,[1 2 150],{'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]},{'full',[0 3]})
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\edited\';file='N7_PreStim_9_8_16';dsf=5;thresh=2.5;onset=5;offset=17000;minInt=5;foi=[1 2 150];eoi={'app',[0 3];'binge',[0 3];'rest',[0 3]};comp=[3];filter='y';overlap=0.5;cycles=3;bands = {'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};

% spectcompbase(sdir,file1,filter,dsf,thresh,onset,offset,minInt,foi,bands,cycles,ftimwin,eventInfo,overlap,comp)

%% Load file
tic
cd(sdir);
disp('Loading file...')
% Load file
load(file);
chans = size(LFPTs.data,1);
toc
%% Filter
if filter == 'y'
    tic
    disp('Applying 60 Hz filter...')
    % Create first order Chebychev stopband filter from 59 to 61 Hz
    [b,a] = cheby1(4,0.5,[59 61]*2/adfreq,'stop');
    % Plot filter
    %fvtool(b2,a2,'Fs',adfreq);
    %figure; plot(LFPTs.data(1,1:20000));
    % Setup filtered LFPT structure
    LFPTsFilt = LFPTs;
    % Apply filter to all chans
    for i = 1:chans
        LFPTsFilt.data(i,:) = filtfilt(b,a,LFPTsFilt.data(i,:));
    end
    % Overwrite LFPTs with LFPTsFilt
    LFPTs = LFPTsFilt;
    toc
else
    disp('Skipping filter...')
    toc
end
%% Downsample data
tic
if dsf > 1
    disp('Downsampling data with dwnSample.m')
    [LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq,chans);
else
    disp('Skipping downsampling...')
end
toc
%% Threshold data
tic
disp('Thresholding data with threshFilt...')
% Got rid of nNaN and indSkp portion
[LFPTs] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq,dsf,chans);
%[LFPTs,nNaN,indSkp] = threshFilt(LFPTs,thresh,onset,offset,minInt,NaNcutoff,adfreq,dsf,chans);
toc
%% Trialize around events
tic
disp('Trializing data with trialize.m')
[eventTs,clnTrls,clnEvents,trls] = trialize(eoi,eventTs,LFPTs,adfreq,minInt,dsf,chans)
toc
%% Calculate power spectra and plot 
tic
disp('Calculating power spectra and plotting average total power with powerComp.m')
if length(comp) == 1
    [psdTrls,relPower,powerPlots] = powerComp(trls,adfreq,eventLabel,chans,comp,bands);
end
if length(comp) == 2
    [psdTrls,relPower,powerPlots,varargout] = powerComp(trls,adfreq,eventLabel,chans,comp,bands);
    stdPower = varargout; clear varargout;
end
toc