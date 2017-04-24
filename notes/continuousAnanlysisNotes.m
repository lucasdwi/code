n = size(LFPTs.data(1,:),2); % 5 secs; 5*400
data = LFPTs.data(1,:);
gamma = 3;
beta = 5;
fs = adfreq;
upperFreq = 100;
upperFreq = upperFreq*(2*pi/fs);
lowerFreq = 1;
lowerFreq = lowerFreq*(2*pi/fs);
f = morsespace(gamma,beta,upperFreq,lowerFreq,4);
psi = morsewave(n,gamma,beta,f);
%%
w = wavetrans(data',psi);
%%
amp = abs(w).^2;
pha = angle(w);
%%
xy = w(:,:,1).*conj(w(:,:,2));
xx = w(:,:,1).*conj(w(:,:,1));
yy = w(:,:,2).*conj(w(:,:,2));
test = (abs(xy).^2)./((abs(xx).^2).*(abs(yy).^2));
%%
autoW = w.^2;
for ii = 1:4
    for jj = 1:4
        crossW(ii,jj,:,:) = w(:,:,ii).*conj(w(:,:,jj));
    end
end
for ii = 1:4
    for jj = 1:4
        coh(ii,jj,:,:) = squeeze(abs(crossW(ii,jj,:,:)).^2./(crossW(ii,ii,:,:).*crossW(jj,jj,:,:)));
    end
end
%% Continuous single-window
% Lowest frequency of interest in Hz
lowFreq = 1;
% Number of cycles per window
minCycle = 5;
% Sampling frequency
fs = adfreq;
% Minimum number of samples for 5 cycles of lowest frequency
samps = (minCycle/lowFreq)*adfreq;
% Next power of 2 to minimum samples
hammWin = 2^nextpow2(samps);
%% Spectrogram
data = LFPTs.data(1,:);

[s,f,t,p] = spectrogram(data,hamming(hammWin),hammWin/2,1:1:100,adfreq);
%% De-noising
thresh = 2.5;
nanDat = LFPTs.data;
nanDat(abs(nanDat)>thresh) = NaN;
nanProp = mean(nanDat,1);
nanDat(:,isnan(nanProp)) = NaN;
%%
start = 1;
c = 1;
while c < floor(length(data)/(hammWin/2))
   avg(c) = mean(nanDat(1,start:start+hammWin));
   start = start + hammWin/2;
   c = c + 1;
end
%% Propogate NaNs to P
pNaN = p;
pNaN(:,isnan(avg)) = NaN;
% Propogate 60 Hz NaNs
pNaN(60,:) = NaN;
%%
figure
h = imagesc(t,f,10*log10(pNaN));
set(h,'alphadata',~isnan(pNaN));
set(gca,'YDir','normal')
