load('paper1data.mat')
coreInd = [11:20,46:50]; % 11:20 pow; 46:50 coh;
shellInd = [1:10,21:25]; % 1:10 pow; 21:25 coh;
powInd = 1:20;
cohInd = 21:50;
thetInd = 1:5:50;
alphInd = 2:5:50;
betInd = 3:5:50;
lgamInd = 4:5:50;
hgamInd = 5:5:50;
% Set up index cell array
inds = {coreInd;shellInd;powInd;cohInd;thetInd;alphInd;betInd;lgamInd;hgamInd};
% Set up index name cell array
indName = {'core','shell','pow','coh','thet','alph','bet','lgam','hgam'};
% Set up cfg for lassoNet
cfg = lassoNetCfg([],'n','y','n',100); 
% Cycle through groups of indices, run lassoNet, and save output
for ii = 1:size(inds,1)
   [allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x(:,inds{ii}),y,'binomial','class',(0:0.01:1),4,1,cfg); 
   save(['C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\subsetLasso\',indName{ii},'.mat'],'allAlpha','allLambda','allBeta','cvFitsArray','accArray','hist')
end
%% Comparisons
files = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\subsetLasso\subsetEnet\'},{'.mat'});

errAvg = [];
for fi = 1:size(files{1},1)
    load(files{1}(fi).name,'allLambda','allBeta');%,'allAlpha')
    beta{fi} = allBeta{1}.betas;
    err{fi} = allLambda{1}.allErr;
    errAvg = [errAvg;allLambda{1}.allErrAvg];
end
%% Core Model: Real-Core-Shell-Rand ECDF
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cvMeanBWHist.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\enetSubsetError.mat')

figure
hold on
ecdf(1-err{4}(:,3))
ecdf(1-err{8}(:,3))
ecdf(1-randData.allErr(:,3))
ecdf(1-realData.allErr(:,3))
title('Core: Subsets ENET')
xlabel('Accuracy')
ylabel('Cumulative Density')
%% Compare Core (4) and Shell (8) features
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cvMeanBWHist.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\subsetError.mat')
%%
subtitles = {'Shell All','Shell Strict','Core All','Core Strict'};
figure
for s = 1:4
   subplot(2,2,s)
   hold on
   ecdf(1-err{4}(:,s))
   ecdf(1-err{8}(:,s))
   ecdf(1-randData.allErr(:,s))
   ecdf(1-realData.allErr(:,s))
   title(subtitles{1,s})
   if s == 1
       legend('Core','Shell','Rand','All')
   end
   xlabel('Accuracy')
   %xlim([.2 1])
end
%% Compare Coh (3) and Power (7) features
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cvMeanBWHist.mat')
subtitles = {'Shell All','Shell Strict','Core All','Core Strict'};
figure
for s = 1:4
   subplot(2,2,s)
   hold on
   ecdf(1-err{3}(:,s))
   ecdf(1-err{7}(:,s))
   ecdf(1-randData.allErr(:,s))
   ecdf(1-realData.allErr(:,s))
   title(subtitles{1,s})
   if s == 1
       legend('Coh','Power','Rand','All')
   end
   xlabel('Accuracy')
   xlim([.2 1])
end
%%
pairs = nchoosek(1:6,2);
ksdat = {realData.allErr;randData.allErr;err{4};err{8};err{3};err{7}};
for ii = 1:size(pairs,1)
    for m = 1:4
        [h(ii,m),p(ii,m)] = kstest2(ksdat{pairs(ii,1)}(:,m),ksdat{pairs(ii,2)}(:,m));
    end
end
pAdj = p.*(numel(p));
%%
for m = 1:4
    for ii = 1:15
        pMat{m}(pairs(ii,1),pairs(ii,2)) = pAdj(ii,m); 
    end
    pMat{m} = pMat{m}';
end

%% Compare Core (4) and Shell (8) Betas
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cvMeanBWHist.mat')
coreSurv = sq(mean(beta{1,4}~=0,1));
coreSurvCut = coreSurv;
coreSurvCut(coreSurvCut<=0.4) = 0;
shellSurv = sq(mean(beta{1,8}~=0,1));
shellSurvCut = shellSurv;
shellSurvCut(shellSurvCut<=0.4) = 0;
%% Run KS tests
for ii = 1:4
    % Core
    [h(1,ii),p(1,ii)] = kstest2(err{4}(:,ii),realData.allErr(:,ii));
    % Shell
    [h(2,ii),p(2,ii)] = kstest2(err{8}(:,ii),realData.allErr(:,ii));
end
%% Compare addition of noise to data
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\subsetLasso\4var.mat','allLambda')
% Color plot of averages
errImage = [allLambda{1}.allErrAvg;errAvg];
figure
imagesc(errImage)
% Plot of cdfs
errDist = cat(3,allLambda{1}.allErr,err
for ii = 1:size(err,2)
    
end
