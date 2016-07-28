%% Open and plot data
name = 'I12BaseSep17.pl2';

PL2tomvdm

for ii = 1:size(LFPTs.data,1)
    figure; plot(LFPTs.data(ii,:))
end