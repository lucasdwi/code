function runSubsetLassoNet(ii)
addpath(genpath('/ihome/ldwiel/code'))
load('/ihome/ldwiel/data/indGroups.mat')
% Set up cfg for lassoNet
cfg = lassoNetCfg([],'n','y','n',100); 
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x(:,inds{ii}),y,'binomial','class',(0:0.01:1),4,1,cfg); 
save(['C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\subsetLasso\',indName{ii},'.mat'],'allAlpha','allLambda','allBeta','cvFitsArray','accArray','hist')