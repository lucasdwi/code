%% Angela Data
%% Reorder channels from Nick's animals (N<20)to be the same as Angela's (N>=20)
LFPTs.data = flipud(LFPTs.data);
LFPTs.label = fliplr(LFPTs.label);
%%
% NAcc (probably) not better than chance at predicting response to stim
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\angelaPFCNAcc.mat')
cfg = lassoNetCfg([],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(NAccX,NAccY,'binomial','class',1,4,1,cfg);
masterRealAcc = 1-allLambda{1,1}.allErr;
% Randomize
cfg = lassoNetCfg([],'y','y','n',100,'1se');
[randAlpha,randLambda,randBeta,randCvFitsArray,randAccArray,randHist] = lassoNet(NAccX,NAccY,'binomial','class',1,4,100,cfg);
masterRandErr = [];
for ii = 1:size(randLambda,2)
    masterRandErr = [masterRandErr;randLambda{1,ii}.allErr];
end
masterRandAcc = 1-masterRandErr;
%%
save('C:\Users\Lucas\Desktop\GreenLab\data\angela\naccData.mat','accArray','allAlpha','allBeta','allLambda','hist','cvFitsArray','masterRandAcc','masterRealAcc','randAccArray','randAlpha','randCvFitsArray','randHist','randLambda')
%% PFC
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\angelaPFCNAcc.mat')
cfg = lassoNetCfg([],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(PFCX,PFCY,'binomial','class',1,4,1,cfg);
masterRealAcc = 1-allLambda{1,1}.allErr;
% Randomize
cfg = lassoNetCfg([],'y','y','n',100,'1se');
[randAlpha,randLambda,randBeta,randCvFitsArray,randAccArray,randHist] = lassoNet(PFCX,PFCY,'binomial','class',1,4,100,cfg);
masterRandErr = [];
for ii = 1:size(randLambda,2)
    masterRandErr = [masterRandErr;randLambda{1,ii}.allErr];
end
masterRandAcc = 1-masterRandErr;
%%
save('C:\Users\Lucas\Desktop\GreenLab\data\angela\pfcData.mat','accArray','allAlpha','allBeta','allLambda','hist','cvFitsArray','masterRandAcc','masterRealAcc','randAccArray','randAlpha','randCvFitsArray','randHist','randLambda')
%%
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\pfcData.mat')
realMean = mean(masterRealAcc);
randMean = mean(masterRandAcc);
figure
ecdf(masterRealAcc)
hold on
ecdf(masterRandAcc)
xlabel('Accuracy')
ylabel('Cumulative Density')
title('Stimulation Response Prediction: PFC')
legend({['Real: ',num2str(round((realMean*100),2)),'%'],['Rand: ',num2str(round((randMean*100),2)),'%']},'location','northwest')
%%
bInd = logicFind(0.2,allBeta{1,1}.survBeta,'>');
load('varNames.mat')
vars = names(bInd);
%%
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\naccData.mat')
realMean = mean(masterRealAcc);
randMean = mean(masterRandAcc);
figure
ecdf(masterRealAcc)
hold on
ecdf(masterRandAcc)
xlabel('Accuracy')
ylabel('Cumulative Density')
title('Stimulation Response Prediction: NAcc')
legend({['Real: ',num2str(round((realMean*100),2)),'%'],['Rand: ',num2str(round((randMean*100),2)),'%']},'location','northwest')
bInd = logicFind(0.2,allBeta{1,1}.survBeta,'>');
load('varNames.mat')
vars = names(bInd);