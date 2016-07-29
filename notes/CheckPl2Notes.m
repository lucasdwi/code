%% Open and plot data
name = 'I1BaseSep14.pl2';

PL2tomvdm;

for ii = 1:size(LFPTs.data,1)
    figure; plot(LFPTs.data(ii,:))
end