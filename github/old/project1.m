%% Binge Data Prep
% Get data structure
hdr = ft_read_header('C:\Users\Lucas\Desktop\PEMM 101\Doucette Lab\Data\H13FoodDepSep7Marked.nex')
% Give clean binge intervals in seconds
% Had to import intervals from Excel spreadsheet used to score intervals on
% how clean the signals were
binge = [22.96643	35.737
37.043	40.217
41.676	52.608
64.76603	97.51
472.5953	481.075
485.734	557.645
558.924	585.7275
2407.501	2419.737
2421.657	2500.582
];
% Convert binge intervals from secons to samples
bingesamp = round(binge * hdr.Fs); 

%% Extract Data from Binges
% Extract Full Relevant Timeserieses
% Only looking at one channel in this version
% N.B.: This data file is not very good due to a hugh notch filter being
% used in preprocessing; regardless, it can be used in this analysis
FP05 = ft_read_data('C:\Users\Lucas\Desktop\PEMM 101\Doucette Lab\Data\H13FoodDepSep7Marked.nex','chanindx',8);
%FP07 = ft_read_data('C:\Users\Lucas\Desktop\PEMM 101\Doucette Lab\Data\H13FoodDepSep7Marked.nex','chanindx',10);
%FP09 = ft_read_data('C:\Users\Lucas\Desktop\PEMM 101\Doucette Lab\Data\H13FoodDepSep7Marked.nex','chanindx',12);
%FP11 = ft_read_data('C:\Users\Lucas\Desktop\PEMM 101\Doucette Lab\Data\H13FoodDepSep7Marked.nex','chanindx',14);
%% Subselect Binge Data
% Should automate
FP05binge = FP05(bingesamp(1,1):bingesamp(1,2));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(2,1):bingesamp(2,2))));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(3,1):bingesamp(3,2))));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(4,1):bingesamp(4,2))));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(5,1):bingesamp(5,2))));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(6,1):bingesamp(6,2))));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(7,1):bingesamp(7,2))));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(8,1):bingesamp(8,2))));
FP05binge = horzcat(FP05binge,(FP05(bingesamp(9,1):bingesamp(9,2))));
%% FFT Plot
np = length(FP05binge);
% Standard periodogram
[PxxP,FP] = periodogram(FP05binge,hamming(np),np,hdr.Fs);
% Welch Spectrum
[PxxW,FW] = pwelch(FP05binge,rectwin(np),np/2,np,hdr.Fs);
%% Plot Both with limited x axis
% Periodogram
subplot(2,2,1)
plot(FP,10*log10(PxxP)); xlim([0 150]);
title('Periodogram')
xlabel('Frequency Hz')
ylabel('Power dB')
% Welch
subplot(2,2,2)
plot(FW,10*log10(PxxW)); xlim([0 150]);
title('Welch')
xlabel('Frequency Hz')
ylabel('Power dB')
%% Smoothing/Averaging
% Lower the frequency resolution to 1Hz steps by averaging across 258
% samples
% N.B.: the Nyquist frequency of the data sampled at 2000 Fs is 1000
% leading to the data being analyzed up to 1 kHz while only the first 150
% Hz are of interest. 

%% Periodogram
chunk = 258;
i = 1;
beg = 1;
PxxPavg = [];
for i = 1:999
    PxxPavg = horzcat(PxxPavg,mean(PxxP(beg:chunk*i)));
    i = i+1;
    beg = beg+chunk;
end
% Downsample frequency axis
FPdwn = FP(1:chunk:length(FP)); 
% takes every 258th element, but this creates a 1000 element vector which
% will not be able to be used with the 999 element vector of PxxPavg. 
% We can take of this when plotting since we are only interested in the
% first 250 Hz anyway.
subplot(2,2,3)
plot(FPdwn(1:150),10*log10(PxxPavg(1:150))); xlim([0 150]);
title('Averaged Periodogram')
xlabel('Frequency Hz')
ylabel('Power dB')
%% Welch
chunk = 258;
i = 1;
beg = 1;
PxxWavg = [];
for i = 1:999
    PxxWavg = horzcat(PxxWavg,mean(PxxW(beg:chunk*i)));
    i = i+1;
    beg = beg+chunk;
end

FWdwn = FW(1:chunk:length(FW));

subplot(2,2,4)
plot(FWdwn(1:150),10*log10(PxxWavg(1:150))); xlim([0 150]);
title('Averaged Welch')
xlabel('Frequency Hz')
ylabel('Power dB')