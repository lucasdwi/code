for ii = 1:size(T1,2)
    binge{1,ii} = T1{1,ii}(:,2:end);
    rest{1,ii} = T2{1,ii}(:,2:end);
    bingeRest{ii} = rest{1,ii} - binge{1,ii};
end
%% Convert grams to kcal
for ii = 1:4
    if ii == 4
        % Multiply grams by kcal/gram in house chow
        bingeCal{ii} = T1{1,ii}(:,1).*3.1;
    else
        % Multiply grams by kcal/gram in sweet-fat food
        bingeCal{ii} = T1{1,ii}(:,1).*4.6;
    end
end
%% Model 1: Baseline Rest vs. Baseline Binge size
cfg = lassoNetCfg([],'n','y','n',100,'min');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(rest{1,1},bingeCal{1,1},'gaussian','mse',1,5,1,cfg);
bInd = logicFind(0.25,allBeta{1,1}.survBeta,'>=');
md1 = fitlm(rest{1,1}(:,bInd),bingeCal{1,1})
%% Model 2: Baseline B-R vs. Baseline Binge Size
cfg = lassoNetCfg([],'n','y','n',100,'min');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(bingeRest{1,1},bingeCal{1,1},'gaussian','mse',1,5,1,cfg);
bInd = logicFind(0.25,allBeta{1,1}.survBeta,'>=');
md2 = fitlm(bingeRest{1,1}(:,bInd),bingeCal{1,1})
%% Model 3: Base->Dep vs. BaseBinge -> Dep24Binge
baseTemp = bingeRest{1,1}(baseDep24(1,:),:);
dep24Temp = bingeRest{1,2}(baseDep24(2,:),:);
xTemp = dep24Temp-baseTemp;
baseBinge = bingeCal{1,1}(baseDep24(1,:),:);
dep24Binge = bingeCal{1,2}(baseDep24(2,:),:);
yTemp = dep24Binge-baseBinge;
cfg = lassoNetCfg([],'n','y','n',100,'min');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(xTemp,yTemp,'gaussian','mse',0:0.1:1,4,1,cfg);
bInd = logicFind(0.25,allBeta{1,1}.survBeta,'>=');
md3 = fitlm(xTemp(:,bInd),yTemp)
%% Model 4: Baseline B-R vs. BaseBinge -> Dep24Binge
xTemp = bingeRest{1,1}(baseDep24(1,:),:);
baseBinge = bingeCal{1,1}(baseDep24(1,:),:);
dep24Binge = bingeCal{1,2}(baseDep24(2,:),:);
yTemp = dep24Binge-baseBinge;
cfg = lassoNetCfg([],'n','y','n',100,'min');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(xTemp,yTemp,'gaussian','mse',1,4,1,cfg);
bInd = logicFind(0.25,allBeta{1,1}.survBeta,'>=');
md4 = fitlm(xTemp(:,bInd),yTemp)
%% Model 5:
% Difference in binge change due to palatability
% Dep24BingeSize - BaseBingeSize
baseDepBinge = bingeCal{1,2}(baseDep24Chow(2,:),:)-bingeCal{1,1}(baseDep24Chow(1,:),:);
% ChowBingeSize - BaseBingeSize
baseChowBinge = bingeCal{1,4}(baseDep24Chow(3,:),:)-bingeCal{1,1}(baseDep24Chow(1,:),:);
% Difference
yTemp = baseDepBinge - baseChowBinge;
% Base B-R - Dep B-R
baseDepEphys = bingeRest{1,1}(baseDep24Chow(1,:),:)-bingeRest{1,2}(baseDep24Chow(2,:),:);
baseChowEphys = bingeRest{1,1}(baseDep24Chow(1,:),:)-bingeRest{1,4}(baseDep24Chow(3,:),:);
xTemp = baseDepEphys - baseChowEphys;
% Lassonet
cfg = lassoNetCfg([],'n','y','n',100,'min');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(xTemp,yTemp,'gaussian','mse',1,3,1,cfg);
bInd = logicFind(0.05,allBeta{1,1}.survBeta,'>=');
md5 = fitlm(xTemp(:,bInd),yTemp)
