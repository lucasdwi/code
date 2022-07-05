figure
subplot(1,2,1)
errorbar(1:7,[0.1,0.5,1.3,1.9,2.4,2.6,3],[0.1,0.3,0.2,0.25,0.4,0.3,0.2],'r')
hold on
errorbar(8:12,[2.9,2.7,2.3,1.4,0.7],[0.35,0.3,0.35,0.25,0.2],'r')
errorbar(1:7,[-0.2,0.3,0.7,1.4,2,2.5,2.3],[0.1,0.15,0.1,0.2,0.4,0.35,0.25],'b')
errorbar(8:12,[1.5,0.9,0.5,0.2,0.3],[0.3,0.15,0.2,0.1,0.05],'b')
errorbar(1:7,[-0.1,0.2,0.15,0.3,0.2,-0.1,0.05],[0.1,0.15,0.1,0.2,0.4,0.35,0.25],'m')
errorbar(8:12,[-0.2,0.12,-0.13,-0.02,0.04],[0.3,0.15,0.2,0.1,0.05],'m')
set(gca,'xtick',1:12,'xticklabel',[1:7,1:5])
ylabel('change in AUC (z-score)')
xlabel('day')
legend({'treatment a','treatment b','treatment c'})
box off

subplot(1,2,2)
hold on
plot(1:7,[0.1,0.18,0.47,0.6,0.65,0.7,0.7],'-or')
plot(8:12,[0.68,0.65,0.5,0.2,0.1],'-or')
plot(1:7,[0.12,0.2,0.49,0.55,0.6,0.65,0.66],'-ob')
plot(8:12,[0.3,0.2,0.1,0.09,0.05],'-ob')
plot(1:7,[0,0,0,0,0,0,0],'-om')
plot(8:12,[0,0,0,0,0],'-om')
set(gca,'xtick',1:12,'xticklabel',[1:7,1:5])
ylabel('% treated')
xlabel('day')
legend({'treatment a','treatment b','treatment c'})
box off
%%
% load('lsd-baseVsal-base_zscore_stimImag_all-216feat.mat')
allBaseZ = zscore(cat(1,allData{1}{:,1}));
allLSD = cat(1,lsdData{1}{[2:7],3});
allSal = cat(1,allData{1}{:,3});
lsdZ = (allLSD-mean(cat(1,allData{1}{:,1}),1))./std(cat(1,allData{1}{:,1}),[],1);
salZ = (allSal-mean(cat(1,allData{1}{:,1}),1))./std(cat(1,allData{1}{:,1}),[],1);
x = 1;
[bN,bE] = histcounts(allBaseZ(:,x),30);
bN = bN./sum(bN);
bE = bE(2:end)-(bE(2)-bE(1))/2;
% histogram(lsdZ(:,x),30)
[sN,sE] = histcounts(salZ(:,x),30);
sE = sE(2:end)-(sE(2)-sE(1))/2;
sN = sN./sum(sN);
figure
hold on
plot(bE,smooth(bN))
plot(sE,smooth(sN))
xlabel('z-score from base')
ylabel('% of samples')
legend({'base','stim'})