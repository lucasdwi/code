%%
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\H10BaseOct15_binge_vs_rest.mat','trls')
%% Generate data
n = 10000; x = randn(n+1,1); data = [x(2:n+1),x(1:n)];
data = data(1:2000,:);
segleng = 128; epleng = 2000;
%%
% Calculate cross-spectral matrix for all channel pairs
[tLen nChan] = size(data);
for ii = 1:nChan
    for k = 1:nChan
        [cs2(ii,k,:),f] = cpsd(data(:,ii),data(:,k),hanning(segleng),segleng/2,400);%,1:200,400);
    end
end
% Calculate complex coherencey
nF = size(cs,3);
for fi = 1:nF
    cc(:,:,fi) = cs(:,:,fi)./sqrt(diag(cs(:,:,fi))*diag(cs(:,:,fi))');
end
% Calculate phase-slope index
df = 1;
psi = imag(sum(conj(cc(:,:,1:end-df)).*cc(:,:,1+df:end),3));
%%
freqbins = [];
[psi2,~,~,~]=data2psi(data,segleng,epleng,freqbins);
%% Concatenate all trials together and transpose
catBingeTrl = [];
for ti = 1:size(trls{1,1}.trial,2)
    catBingeTrl = [catBingeTrl;trls{1,1}.trial{1,ti}'];
end
catRestTrl = [];
for ti = 1:size(trls{1,2}.trial,2)
    catRestTrl = [catRestTrl;trls{1,2}.trial{1,ti}'];
end

%% Run PSI
% Generate freq axis for binning
f = 400*(0:1000)/2000;
ends = nearest_idx3([1:100],f);
starts = [2;ends(1:end-1)+1];
for ii = 2:size(ends,1)
    freqs(ii,:) = starts(ii):ends(ii);
end
%%
freqs = [];
freqs(:,1) = 1:2:500;
freqs(:,2) = 2:2:500;
%%
fStart = nearest_idx3(45,f);
fStop = nearest_idx3(65,f);
[psi,stdpsi,psisum,stdpsisum] = data2psi(catBingeTrl,2000,2000,freqs);
normPsiBinge = psi./(stdpsi+eps);
[psi,stdpsi,psisum,stdpsisum] = data2psi(catRestTrl,2000,2000,freqs);
normPsiRest = psi./(stdpsi+eps);
%% Freq Plot
figure
hold on
plot(1:250,squeeze(normPsiBinge(1,3,:)))
plot(1:250,squeeze(normPsiRest(1,3,:)))
%% Plot
% Get max and min for colormap
cMax = max(max(normPsiBinge))
figure
subplot(1,2,1)
pcolor(padarray(normPsiBinge,[1 1],0,'post'))
set(gca,'XTick',1.5:1:4.5,'XTickLabel',1:4,'YTick',1.5:1:4.5,'YTickLabel',1:4)
subplot(1,2,2)
pcolor(padarray(normPsiRest,[1 1],0,'post'))
set(gca,'XTick',1.5:1:4.5,'XTickLabel',1:4,'YTick',1.5:1:4.5,'YTickLabel',1:4)
colormap('viridis')

