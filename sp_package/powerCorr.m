function [r1,r2] = powerCorr(psdTrls,bands)
%% Calculates power correlations within frequency bands across channels
% INPUTS:
% psdTrls = 
% bands = 
%__________________________________________________________________________
% OUTPUTS:
% r1 = 
% r2 = 
%__________________________________________________________________________
%% LLD and MAC 2017
%% Preallocate different sizes
chans = size(psdTrls.event1.Overall,1);
nTrls = [size(psdTrls.event1.Pow,2),size(psdTrls.event2.Pow,2)];
nBands = size(bands,1);
% Grab frequency band indices using frequency vector from powerComp.m
bInd = bandIndices(bands,psdTrls.F);
% Event 1 band power normalization
event1Pow = zeros(nBands,chans,nTrls(1,1));
total1Pow = zeros(nBands,chans,nTrls(1,1));
for ii = 1:size(psdTrls.event1.Pow(2,:),2)
    event1Pow(:,:,ii) = psdTrls.event1.Pow{2,ii};
    total1Pow(:,:,ii) = repmat(trapz(psdTrls.event1.Pow{1,ii}(:,bInd(1,1):bInd(end,2)),2)',6,1,1);
end
norm1Pow = event1Pow./total1Pow;
% Event 2 band power normalization
event2Pow = zeros(nBands,chans,nTrls(1,2));
total2Pow = zeros(nBands,chans,nTrls(1,2));
for ii = 1:size(psdTrls.event2.Pow(2,:),2)
    event2Pow(:,:,ii) = psdTrls.event2.Pow{2,ii};
    total2Pow(:,:,ii) = repmat(trapz(psdTrls.event2.Pow{1,ii}(:,bInd(1,1):bInd(end,2)),2)',6,1,1);
end
norm2Pow = event2Pow./total2Pow;
%
r1 = zeros(chans,chans,nBands); p1 = zeros(size(r1));
r2 = zeros(chans,chans,nBands); p2 = zeros(size(r2));
for iB = 1:nBands
    [r1(:,:,iB),p1(:,:,iB)] = corrcoef(squeeze(norm1Pow(iB,:,:))');
    [r2(:,:,iB),p2(:,:,iB)] = corrcoef(squeeze(norm2Pow(iB,:,:))');
end
