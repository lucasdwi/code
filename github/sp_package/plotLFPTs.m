function plotLFPTs(LFPTs,adfreq)
%% Converts tvec into secs and plots
plotTvec = 0:length(LFPTs.tvec)/adfreq/60;
figure;
subplot(2,2,1)
plot(plotTvec,LFPTs.data(1,:))
subplot(2,2,2)
plot(plotTvec,LFPTs.data(2,:))
subplot(2,2,3)
plot(plotTvec,LFPTs.data(3,:))
subplot(2,2,4)
plot(plotTvec,LFPTs.data(4,:))