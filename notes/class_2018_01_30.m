%% FILTER
load('I1FoodDep24_2015-11-13.mat')
data = LFPTs; 
%%
[data.data,data.tvec] = sigGen(2000,10,[4,60],[1,1.5],0);
noise = normrnd(0,5,[1,20001]);
data.data = data.data + noise;
%%
figure
subplot(4,1,1)
plot(data.tvec(1:20001),data.data(1,1:20001))
hold on
% plot(data.tvec,sin(2*pi*4*data.tvec))
% plot(data.tvec,sin(2*pi*60*data.tvec))
xlim([0 2])
xlabel('Time (sec)')
ylabel('mV')
%%
[x,f] = pwelch(data.data(1,4000*2000:4200*2000),2048,[],1:100,2000);
subplot(4,1,2)
plot(f,10*log10(x))
xlabel('Frequency (Hz)')
ylabel('dB')
%%
filtData = filter60(data.data(1,:),2000,0);
subplot(4,1,3)
plot(data.tvec,filtData);
hold on
% plot(data.tvec,sin(2*pi*4*data.tvec)*10)
% plot(data.tvec,sin(2*pi*60*data.tvec)*10)
xlim([0 2])
xlabel('Time (sec)')
ylabel('mV')
%%
[x,f] = pwelch(filtData(1,1:20001),2048,[],1:100,2000);
subplot(4,1,4)
plot(f,10*log10(x))
xlabel('Frequency (Hz)')
ylabel('dB')
%% DOWNSAMPLE
%%
figure
subplot(2,1,1)
plot(data.tvec,filtData)
xlim([0 0.25])
xlabel('Time (sec)')
ylabel('mV')
%%
LFPTs.data = filtData;
[LFPTs,adfreq] = dwnSample(LFPTs,5,adfreq);
subplot(2,1,2)
plot(LFPTs.tvec,LFPTs.data)
xlim([0 0.25])
xlabel('Time (sec)')
ylabel('mV')
%%
load('H10Base_2015-09-24.mat')