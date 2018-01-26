%% Fake data
% Creates data with two 'behaviors' with separate frequencies; b1 should
% have 12 complete 5 second trials with a 3 Hz component; b2 should have 3
% complete 5 second trials at 11 Hz; neither 63 Hz (noise) nor 74 Hz (not
% enough data) should be present
fs = 1000;
int = 5;
% Channel 1
[s1,~] = sigGen(fs,int*10,3,1,0);
[s2,~] = sigGen(fs,50.017,63,1,0);
noise = [0:1/4:4,4*ones(1,5000),-4*ones(1,5000),-4:1/10000:0];
s2 = s2.*noise;
[s3,~] = sigGen(fs,int*2.5,3,1,0);
[s4,~] = sigGen(fs,int*3,11,1,0);
[s5,~] = sigGen(fs,int*0.8,74,1,0);
x1 = cat(2,s1,s2,s3,s4,s2,s5);
% plot(x1)
% Channel 2, copy of channel one but with slight phase shift and more clean
% data
[s1,~] = sigGen(fs,int*11.6,3,1,10);
[s2,~] = sigGen(fs,42.017,63,1,10);
noise = [0:1/4:4,4*ones(1,1000),-4*ones(1,1000),-4:1/10000:0];
s2 = s2.*noise;
[s3,~] = sigGen(fs,int*2.5,3,1,10);
[s4,~] = sigGen(fs,int*4.6,11,1,10);
[s5,~] = sigGen(fs,int*0.8,74,1,10);
x2 = cat(2,s1,s2,s3,s4,s2,s5);
% hold on
% plot(x2)

eventTs.label = [{'b1 (Start)','b1 (End)','b2 (Start)','b2 (End)'}];
eventTs.t = {[0;100.019],[50.001;112.52],[112.521,177.539],[135.521,181.54]};

LFPTs.data = cat(1,x1,x2);
LFPTs.label = {'Channel 1','Channel 2'};
LFPTs.tvec = linspace(0,181.54,length(LFPTs.data));
LFPTs.cfg.history.mfun = [];
adfreq = fs;
%%

%%