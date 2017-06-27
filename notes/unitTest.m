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
time = round(rand(1)*3e3);
fs = 1000;
[x,tvec] = sigGen(fs,time,60,1,[]);
LFPTs.data = repmat(x,4,1);
LFPTs.tvec = tvec;
LFPTs.cfg.history.mfun = [];
%% Test DSF 1
[LFPTs,adfreq] = dwnSample(LFPTs,1,fs);
assert(length(LFPTs.data) == length(x)/1)
%% Test DSF 2
[LFPTs,adfreq] = dwnSample(LFPTs,2,fs);
assert(length(LFPTs.data) == length(x)/2)
%% Test DSF 3
[LFPTs,adfreq] = dwnSample(LFPTs,3,fs);
assert(length(LFPTs.data) == length(x)/3)
