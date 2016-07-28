%%
% load file and create pl2 struct
cd('C:\Users\Lucas\Desktop');
fname = name;
pl2 = PL2GetFileIndex(name);

%%
% Display all ad channels in file
PL2Print(pl2.AnalogChannels);
%%
%check for block start time and sample number per block

%bdata = nan(length(pl2.AnalogChannels{1, 1}.InternalArray1)-1,2);
%cblock = PL2ReadFirstDataBlock(fname);
%bdata(1,1) = cblock.AnalogData.NumSamples;
%bdata(1,2) = cblock.AnalogData.Timestamp;
%for b = 2:length(pl2.AnalogChannels{1, 1}.InternalArray1)-1
%    cblock = internalPL2ReadNextDataBlock(cblock);
%    bdata(b,1) = cblock.AnalogData.NumSamples;
%    bdata(b,2) = cblock.AnalogData.Timestamp;
%end
%% add ad channel frequency data
[n, freqs] = plx_adchan_freqs(fname); 
%%
% Load lfp ad data and labels into tsd format struct called LFPTs
    %Remove all empty and wide band ad channels
    lfpchan = nan(length(pl2.AnalogChannels), 1);
for j = 1:length(pl2.AnalogChannels); %create logical of good channels
    lfpchan(j,1) = pl2.AnalogChannels{j,1}.NumValues > 0 & freqs(j,1)<5000;
end
lfpchan = find(lfpchan); %convert logical to indexes of good channels
%%
 if (isempty(lfpchan)) == 0  
    LFPTs = tsd([]);
    
    LFPTs.data = nan(length(lfpchan),pl2.AnalogChannels{lfpchan(1),1}.NumValues);%empty data array
    for i = 1:length(lfpchan)%fill data array from all the lfpchan
        LFPTs.label{1,i} = pl2.AnalogChannels{lfpchan(i),1}.Name;
        temp = PL2Ad(fname, lfpchan(i));
        LFPTs.data(i,:) = temp.Values;
    end
    
    % Add time vector to LFPTs.tvec
    [adfreq, n, ts, fn, ad] = plx_ad(fname, lfpchan(1));
    LFPTs.tvec = nan(1, pl2.AnalogChannels{lfpchan(1),1}.NumValues);
    LFPTs.tvec(1,:) = ts:1/freqs(lfpchan(1)):(pl2.AnalogChannels{lfpchan(1),1}.NumValues-1)/freqs(lfpchan(1))+ts;
    TimeSampEr = PL2StartStopTs(fname, 'stop')- LFPTs.tvec(length(LFPTs.tvec));
 else
 end
    %%
    % Extraction of wideband
    %%
    % Load WB ad data and labels into tsd format struct called WBTs
        %Remove all empty and lfp ad channels
        WBchan = nan(length(pl2.AnalogChannels), 1);
    for j = 1:length(pl2.AnalogChannels); %create logical of good channels
        WBchan(j,1) = pl2.AnalogChannels{j,1}.NumValues > 0 & freqs(j,1)>5000;
    end
    WBchan = find(WBchan); %convert logical to indexes of good channels
%%
if (isempty(WBchan)) == 0 

    WBTs = tsd([]);
    WBTs.data = nan(length(WBchan),pl2.AnalogChannels{WBchan(1),1}.NumValues);%empty data array
    for i = 1:length(WBchan)%fill data array from all the WBchan
        WBTs.label{1,i} = pl2.AnalogChannels{WBchan(i),1}.Name;
        temp = PL2Ad(fname, WBchan(i));
        WBTs.data(i,:) = temp.Values;
    end

% Add time vector to WBTs.tvec
    [adfreq, n, ts, fn, ad] = plx_ad(fname, WBchan(1));
    WBTs.tvec = nan(1, pl2.AnalogChannels{WBchan(1),1}.NumValues);
    WBTs.tvec(1,:) = ts:1/freqs(WBchan(1)):(pl2.AnalogChannels{WBchan(1),1}.NumValues-1)/freqs(WBchan(1))+ts;
    TimeSampEr = PL2StartStopTs(fname, 'stop')- WBTs.tvec(length(WBTs.tvec));   
else
end
%%
%remove nans
%AdTs.data(isnan(AdTs.data))=0; 
%%
% Display all event channels in file
PL2Print(pl2.EventChannels)
%%
% get all event timestamps and labels from PL2 file and put into ts struct
eventTs = ts([]);
for i = 1:length(pl2.EventChannels)
    eventTs.label{1,i} = pl2.EventChannels{i,1}.Name;
    temp = PL2EventTs(fname, i);
    eventTs.t{1,i} = temp.Ts;
end
