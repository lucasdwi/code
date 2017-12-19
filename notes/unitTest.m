% %% Creates signal with 60 cycle line noise and filters it out
% [x,tvec] = sigGen(1000,120,60,1,[]);
% x = x+randi([-1 1],1,length(x)).*randn(1,length(x));
% LFPTs.data = repmat(x,4,1);
% [filtData] = filter60(LFPTs,1000,'n');
% [pxxFilt] = pwelch(filtData',256,128,1:100,1000);
% [pxx] = pwelch(LFPTs.data',256,128,1:100,1000);
% % Plots results
% figure
% subplot(1,2,1)
% plot(tvec,x,'b-')
% hold on
% plot(tvec,filtData(1,:),'r-')
% xlim([0 2])
% legend({'Original','Filtered'},'location','north')
% subplot(1,2,2)
% plot(1:100,10*log10(pxx),'b-');
% hold on
% plot(1:100,10*log10(pxxFilt),'r-');
% Test dwnSample.mat; checks that number of samples is factor (dsf) of original
%% Test dwnSample
time = 3456.789;
fs = 1000;
[x,tvec] = sigGen(fs,time,60,1,[]);
LFPTs.data = repmat(x,4,1);
LFPTs.tvec = tvec;
LFPTs.cfg.history.mfun = [];
dsf = 5;
% Run test
[testData,adfreq] = dwnSample(LFPTs,dsf,fs);
assert(length(testData.data) == length(x)/dsf && length(testData.tvec) == length(testData.data))
assert(adfreq == fs/dsf)
%% Test eventInd
% Set up  eventTs
eventTs.label = {'Int (Start)','Int (End)','Scalar'};
eventTs.t = {1,10,15;21,30,35};
eoi = {'sca',[-1 1];'int',[0 5]};
[eventInds] = eventInd(eventTs,eoi);
assert(isequal(eventInds,[3,3;1,2]))
eoi = {'int (e',[-1 1]};
[eventInds] = eventInd(eventTs,eoi);
assert(isequal(eventInds,[2,2]))
%% Test threshFilt
load('fakeData.mat')
% Run test
thresh = 2; onset = 0.011; offset = 20; minInt = 5;
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq);
nanPosition = logicFind(1,isnan(LFPTs.data(1,:)),'==');
diffNaN = diff(nanPosition);
ind = logicFind(1,diffNaN,'~=');
start = nanPosition([1,ind+1]);
stop = nanPosition([ind,end]);
% Tests starts and stops of NaNed data
assert(isequal(start,[50002,127522]))
assert(isequal(stop,[100006,181540]))
% Tests that no channels were removed
assert(isempty(zeroedChannel))
%% Test trialize
load('fakeData.mat')
thresh = 2; onset = 0.011; offset = 20; minInt = 5;
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq);
% eoi of both 'behaviors'
eoi = {'b1',[0 5]; 'b2', [0 5]}; 
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
assert(size(trls{1,1}.trial,3)==12)
assert(size(trls{1,2}.trial,3)==3)
% all
eoi = {'all',[0 5]}; 
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
assert(size(trls{1,1}.trial,3)==15)
% notbinge
% eoi = {'not',[0 5]}; 
% [clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
% assert(size(trls{1,1}.trial,3)==15)
%% Test powerComp
load('fakeData.mat')
thresh = 2; onset = 0.011; offset = 20; minInt = 5;
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq);
eoi = {'b1',[0 5]; 'b2', [0 5]}; 
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);

chans = size(LFPTs.data,1);
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
filter = 'n';
foi = [1 2 500];
[psdTrls,~,hist] = powerComp(trls,adfreq,chans,bands,filter,foi,eoi,'n');
figure
subplot(2,2,1)
plot(hist.powF,psdTrls{1,1}.Overall)
title('pwelch: b1 3 Hz')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
subplot(2,2,2)
plot(hist.powF,psdTrls{1,2}.Overall)
title('pwelch: b2 11 Hz')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
% Grab peaks of b1 power
[~,inds] = max(psdTrls{1,1}.Overall,[],2);
freq1 = hist.powF(inds);
assert(isequal(freq1,[3,3]))
% Grab peaks of b2 power
[~,inds] = max(psdTrls{1,2}.Overall,[],2);
freq2 = hist.powF(inds);
assert(isequal(freq2,[11,11]))

% Check using fft as well
for iB = 1:size(eoi,1)
    for iT = 1:size(trls{1,iB}.trial,3)
        pow{iB}(:,:,iT) = fft(trls{1,iB}.trial(:,:,iT),5000,2);
    end
    mag(:,:,iB) = mean(abs(pow{iB}/2500),3);
end
mag = mag(:,1:2501,:);
f = linspace(0,500,2501);
fInds = [];
for fi = [3,11,63,74]
    fInds = [fInds,logicFind(fi,f,'==')];
end
% Round to get rid of trailing zeros; check against known amplitudes for
% 3,11,63, and 74 Hz in the two behaviors
assert(isequal(round(mag(:,fInds,:)),cat(3,[1,0,0,0;1,0,0,0],[0,1,0,0;0,1,0,0])))
subplot(2,2,3)
plot(f,mag(:,:,1))
title('fft: b1 3 Hz')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
subplot(2,2,4)
plot(f,mag(:,:,2))
title('fft: b2 11 Hz')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
tightfig;
%% Test cohComp
load('fakeData.mat')
thresh = 2; onset = 0.011; offset = 20; minInt = 5;
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq);
eoi = {'b1',[0 5]; 'b2', [0 5]}; 
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);

chans = size(LFPTs.data,1);
cycles = 5;
ftimwin = [];
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
foi = [1 2 500];
overlap = 0.5;
if ~isempty(ftimwin)
    [~,winSize] = nearestPow2(ftimwin*adfreq);
elseif ~isempty(cycles)
    % Uses lowest frequency (first in bands) to compute winSize
    [~,winSize] = nearestPow2((cycles/bands{1,2}(1))*adfreq);
end
[coh,cohPlots] = cohCompMat(LFPTs,chans,trls,foi,winSize,overlap,adfreq,bands,zeroedChannel,eoi);
