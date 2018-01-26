

winSizes = 0.4:.2:5;
overlaps = .1:.1:0.9;
for ii = 1:length(winSizes)
    [~,winSize] = nearestPow2(adfreq*winSizes(ii));
    window = hamming(winSize);
    for jj = 1:length(overlaps)
        overlapSamp = round(overlaps(jj)*winSize);
        for ti = 1:size(trls{1}.trial,3)
            x = trls{1}.trial(1,:,ti);
            y = trls{1}.trial(2,:,ti);
            [Cxy{ii,jj}(:,ti),f] = mscohere(x,y,window,overlapSamp,foiV,...
                adfreq);
        end
    end
    disp([num2str(ii),' of ',num2str(length(winSizes))])
end
%%
for mi = 1:200
    x1 = randn(1,10001);
    d = designfilt('lowpassiir','passbandfrequency',5,'StopbandFrequency',150,'PassbandRipple',0.5,'StopbandAttenuation',65,'DesignMethod','butter','SampleRate',1000);
    d2 = designfilt('highpassiir','passbandfrequency',25,'StopbandFrequency',10,'PassbandRipple',0.5,'StopbandAttenuation',65,'DesignMethod','butter','SampleRate',1000);
    x2 = filter(d,x1)+normrnd(0,0.1,1,10001);
    [h,w] = freqz(d);
    
    trueCoh(mi,:) = abs(h).^2./(abs(h).^2+0.1);
   
%     wSizes = [64,128,256,512,1024,2048,4096];
    wSizes = 256:4:512;
    overlaps = 0.1:0.1:0.9;
    for ii = 1:size(wSizes,2)
        for jj = 1:size(overlaps,2)
            [Cxy(ii,jj,:,mi),f] = mscohere(x1,x2,hamming(wSizes(ii)),round(overlaps(jj)*wSizes(ii)),1:1:100,1000);
        end
    end
end
%%
f2 = (w*1000)./(2*pi);
figure
c = 1;
for ii = 21:30
    for jj = 1:9
        subplot(10,9,c)
        plot(f,squeeze(mean(Cxy(ii,jj,:,:),4)))
        hold on
        plot(f2(1:103),trueCoh(1,1:103))
        c = c+1;
    end
end