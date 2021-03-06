%% Grab data
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkNot\',{'.mat'},{'pow','coh'},'trl','rel');
% Concat all together to avoid issues of low drinking samples 
catData{1,1} = cat(1,data{1,1}{:,1});
catData{1,2} = cat(1,data{1,1}{:,2});
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drinkNotRel.mat'],'catData')
%% Run evenDataSplit: use 50-50 training, but force train/test into 80/20
[all,~,rnd,~] = evenDataSplit(catData,1000,250,'ADA',100);
% [all,each,rnd,weights] = evenDataSplit(catData,6000,1200,'ADA',20);
%% Save
save('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\concatData1000.mat','all','rnd')
% save('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drink6000.mat','all','rnd')
%% Grab water data
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\waterProcessed\'],{'.mat'},{'pow','coh'},'trl','rel');
waterData{1,1} = cat(1,data{1,1}{:,1});
waterData{1,2} = cat(1,data{1,1}{:,2});
[allWater,~,rndWater,~] = evenDataSplit(waterData,1000,250,'ADA',100);
%% Determine relationship between number of samples per animal, # of recordings, and total # of samples
alcN = cellfun(@(x) size(x,1),alcData(:,1),'UniformOutput',0);
watN = cellfun(@(x) size(x,1),watData(:,1),'UniformOutput',0);
for ii = 1:50
    nAlc(ii) = sum(cell2mat(alcN)>=ii);
    nWat(ii) = sum(cell2mat(watN)>=ii);
    
    nAlcSamp(ii) = ii*nAlc(ii);
    nWatSamp(ii) = ii*nWat(ii);
end
figure
hold on
plot(nAlc,nAlcSamp)
plot(nWat,nWatSamp)
figure
hold on
plot3(1:ii,nAlc,nAlcSamp)
plot3(1:ii,nWat,nWatSamp)
plot3(1:ii,nWat+nAlc,nWatSamp+nAlcSamp)
xlabel('Samp per animal')
ylabel('# Recordings')
zlabel('# Samp')
%% Only pull 10 trials from each animal
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\alcData.mat')
alcData = data;
alcFiles = files;
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\waterData.mat')
waterData = data;
waterFiles = files;
% Get indices of zero drinking and >= 10 drinking
zeroIndAlc = logicFind(0,cellfun(@(x) size(x,1),alcData{1,1}(:,1)),'==');
greaterIndAlc = logicFind(10,cellfun(@(x) size(x,1),alcData{1,1}(:,1)),'>=');

zeroIndWater = logicFind(0,cellfun(@(x) size(x,1),waterData{1,1}(:,1)),'==');
greaterIndWater = logicFind(10,cellfun(@(x) size(x,1),waterData{1,1}(:,1)),'>=');
% Water has 2 pairs and 2 singles with zero; 3 pairs, 7 singles
% with >=10
alcGreaterPairs = {3:6;7:10;11:15};

alcZeroPairs = {1:2;3:6};
%% Manually pull data from equivalent number of animals/recordings
% Grab 10 bins from all greaterWater data
for jj = 1:100
    thisDrink = [];
    thisNot = [];
    for ii = greaterIndWater
        thisDrink = [thisDrink;...
            waterData{1,1}{ii,1}(randi(size(waterData{1,1}{ii,1},1),10,1),:)];
        thisNot = [thisNot;...
            waterData{1,1}{ii,2}(randi(size(waterData{1,1}{ii,2},1),700,1),:)];
    end
    waterCat{1,1} = thisDrink;
    waterCat{1,2} = thisNot;
    [thisAll,~,thisRnd,~] = evenDataSplit(waterCat,1000,250,'ADA',1);
    allWater.trainX(jj) = thisAll.trainX;
    allWater.trainY(jj) = thisAll.trainY;
    allWater.testX(jj) = thisAll.testX;
    allWater.testY(jj) = thisAll.testY;
    
    rndWater.allTrainY(jj) = thisRnd.allTrainY;
    rndWater.allTestY(jj) = thisRnd.allTestY;
    
    
    thisDrink = [];
    thisNot = [];
    for ii = [3,10,28,30,33,34]
        thisDrink = [thisDrink;...
            alcData{1,1}{ii,1}(randi(size(alcData{1,1}{ii,1},1),10,1),:)];
        thisNot = [thisNot;...
            alcData{1,1}{ii,2}(randi(size(alcData{1,1}{ii,2},1),700,1),:)];
    end
    alcCat{1,1} = thisDrink;
    alcCat{1,2} = thisNot;
    [thisAll,~,thisRnd,~] = evenDataSplit(alcCat,1000,250,'ADA',1);
    all.trainX(jj) = thisAll.trainX;
    all.trainY(jj) = thisAll.trainY;
    all.testX(jj) = thisAll.testX;
    all.testY(jj) = thisAll.testY;
    
    rndWater.allTrainY(jj) = thisRnd.allTrainY;
    rndWater.allTestY(jj) = thisRnd.allTestY;
end

%% Build models - no 'each' models due to cases of low/no drinking samples 
% Preallocate
% [x,y,t,xRand,yRand] = deal(cell(1,20));
% [a,aRand] = deal(zeros(1,20));
for n = 1:100
    disp([num2str(n),' of 100'])
    % Set up training data
    trainX = all.trainX{n};
    trainY = all.trainY{n};
    % Set up testing data
    testX = all.testX{n};
    testY = all.testY{n};
    % Build and test model on concat data
    aa.mdl{n} = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(aa.mdl{n},testX);
    [aa.x{n},aa.y{n},~,aa.a(n)] = perfcurve(testY,prob,1,...
        'TVals',0:1/4623:1,'UseNearest',0);
    % Apply to random
    trainY = rnd.allTrainY{n};
    testY = rnd.allTestY{n};
    prob = predict(aa.mdl{n},testX);
    [aa.rndX{n},aa.rndY{n},~,aa.rndA(n)] = perfcurve(testY,...
        prob,1,'Tvals',0:1/4623:1,'UseNearest',0);
end
%% Calculate average
mX = mean(cat(2,aa.x{:}),2);
mY = mean(cat(2,aa.y{:}),2);
mRandX = mean(cat(2,aa.rndX{:}),2);
mRandY = mean(cat(2,aa.rndY{:}),2);
%%
figure
hold on
plot(mX,mY,'-k');
plot(mRandX,mRandY,'--k');
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
ylabel('True Positive Rate')
xlabel('False Positive Rate')
legend({['Real: ',num2str(round(mean(aa.a),2))],...
    ['Permuted: ',num2str(round(mean(aa.rndA),2))]},'Location','se')
title('Drink vs. Other')
%% Water v. alcohol
% Cross compare feature lists
waterFeat = names({'ILL','CA1L','PL','SL','PR','CA1R','ILR','SR'},{'d','t','a','b','lg','hg'});
alcFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
% Find overlap
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterFeat,'==');
    if ~isempty(ind)
        inds(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
inds(43:48) = 157:162;
%% Build using above indices for water (inds)
for n = 1:100
    disp([num2str(n),' of 100'])
    % Set up training data
    alcTrainX = all.trainX{n}{1};
    alcTrainY = all.trainY{n}{1};
    % Set up testing data
    alcTestX = all.testX{n}{1};
    alcTestY = all.testY{n}{1};
    % Build and test alcohol model
    aa.mdl{n} = fitglm(alcTrainX,alcTrainY,'distribution','binomial');
    prob = predict(aa.mdl{n},alcTestX);
    [aa.x{n},aa.y{n},~,aa.a(n)] = perfcurve(alcTestY,prob,1,...
        'TVals',0:1/4623:1,'UseNearest',0);
    % Apply alcohol model to random
    rndTestY = rnd.allTestY{n};
    prob = predict(aa.mdl{n},alcTestX);
    [aa.rndX{n},aa.rndY{n},~,aa.rndA(n)] = perfcurve(...
        rndTestY,prob,1,'Tvals',0:1/4623:1,'UseNearest',0);
    % Apply alcohol model to water
    waterTestX = allWater.testX{n}{1}(:,inds);
    waterTestY = allWater.testY{n}{1};
    prob = predict(aa.mdl{n},waterTestX);
    [aw.x{n},aw.y{n},~,aw.a(n)] = perfcurve(...
        waterTestY,prob,1,'TVals',0:1/4623:1,'UseNearest',0);
    % Build water model
    waterTrainX = allWater.trainX{n}{1}(:,inds);
    waterTrainY = allWater.trainY{n}{1};
    ww.mdl{n} = fitglm(waterTrainX,waterTrainY,'distribution','binomial');
    prob = predict(ww.mdl{n},waterTestX);
    [ww.x{n},ww.y{n},~,ww.a(n)] = perfcurve(waterTestY,prob,1,'TVals',...
        0:1/4623:1,'UseNearest',0);
    % Apply water model to alcohol
    prob = predict(ww.mdl{n},alcTestX);
    [wa.x{n},wa.y{n},~,wa.a(n)] = perfcurve(alcTestY,prob,1,'TVals',...
        0:1/4623:1,'UseNearest',0);
end
%%
figure
hold on
plot(mean(cat(2,aa.x{:}),2),mean(cat(2,aa.y{:}),2),'-b')
plot(mean(cat(2,aw.x{:}),2),mean(cat(2,aw.y{:}),2),'-c')
plot(mean(cat(2,ww.x{:}),2),mean(cat(2,ww.y{:}),2),'-r')
plot(mean(cat(2,wa.x{:}),2),mean(cat(2,wa.y{:}),2),'-m')
plot(mean(cat(2,aa.rndX{:}),2),mean(cat(2,aa.rndY{:}),2),'--k')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
ylabel('True Positive Rate')
xlabel('False Positive Rate')
legend({['A-A: ',num2str(round(mean(aa.a),2)),'\pm',...
    num2str(round(conf(aa.a,0.95),2))],['A-W: ',...
    num2str(round(mean(aw.a),2)),'\pm',...
    num2str(round(conf(aw.a,0.95),2))],['W-W: ',...
    num2str(round(mean(ww.a),2)),'\pm',...
    num2str(round(conf(ww.a,0.95),2))],['W-A: ',...
    num2str(round(mean(wa.a),2)),'\pm',...
    num2str(round(conf(wa.a,0.95),2))],['Rnd: ',...
    num2str(round(mean(aa.rndA),2)),'\pm',...
    num2str(round(conf(aa.rndA,0.95),2))]},'location','se')
title('Alcohol vs. Water')
%%
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'logCmbs_concatData\1.mat'])
monadM = mean(A,1);
monadS = std(A,[],1);
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
[sortMonadM,ord] = sort(monadM','descend');
sortMonadS = monadS(ord)';
monadFeats = nameVect(ord)';
monadTier = tier(A(:,ord));
%% Plot
figure
% shadedErrorBar(1:60,mean(A(:,ord),1),std(A(:,ord),[],1))
scatterErr(1:60,mean(A(:,ord),1),std(A(:,ord),[],1),0)
hold on
for ii = 1:size(monadTier,2)
   plot([monadTier(ii) monadTier(ii)],[0 1],'--k','LineWidth',1) 
end
ylim([0.3 1])
%%
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\logCmbs\2.mat')
dyadM = mean(A,1);
dyadS = std(A,[],1);
[sortDyadM,ord] = sort(dyadM','descend');
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
cmbs = nchoosek(1:60,2);
for ii = 1:size(cmbs,1)
   dyadFeats(ii,1) = nameVect(cmbs(ord(ii),1));
   dyadFeats(ii,2) = nameVect(cmbs(ord(ii),2));
end
sortDyadS = dyadS(ord)';
dyadTier = tier(A(:,ord));
%%
figure
shadedErrorBar(1:1770,mean(A(:,ord),1),std(A(:,ord),[],1))
hold on
for ii = 1:size(dyadTier,2)
   plot([dyadTier(ii) dyadTier(ii)],[0 1],'--k','LineWidth',1) 
end
ylim([0.3 1])
xlim([0 1770])
%%
for ii = 1:20
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\logCmbs\triplets\trip',num2str(ii),'.mat']) 
   a(ii,:) = A;
end
triadM = mean(a,1);
triadS = std(a,[],1);
[sortTriadM,ord] = sort(triadM','descend');
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
cmbs = nchoosek(1:60,3);
for ii = 1:size(cmbs,1)
   triadFeats(ii,1) = nameVect(cmbs(ord(ii),1));
   triadFeats(ii,2) = nameVect(cmbs(ord(ii),2));
   triadFeats(ii,3) = nameVect(cmbs(ord(ii),3));
end
sortTriadS = triadS(ord)';
triadTier = tier(a(:,ord));
%%
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drinkNotModel.mat','concatData')
scatterErr(1:4,[monadM(1),dyadM(1),triadM(1),mean(concatData.auc)],[monadS(1),dyadS(1),triadS(1),std(concatData.auc,[],2)],1)
set(gca,'XTick',1:4,'XTickLabel',{'Monad','Dyad','Triad','Logistic (60)'})
ylabel('AUC')
title('Drink vs. Other: Feature #')
xlim([0.75 4.25])
%%
freqs = {'d','t','a','b','lg','hg'};
for ii = 1:size(freqs,2)
   perc(ii,1) = sum(~cellfun(@isempty,regexp(monadFeats(1:(monadTier(end)-1),:),freqs{ii})))/(monadTier(end)-1);
   perc(ii,2) = sum(any(~cellfun(@isempty,regexp(dyadFeats(1:(dyadTier(end)-1),:),freqs{ii})),2))/(dyadTier(end)-1);
   perc(ii,3) = sum(any(~cellfun(@isempty,regexp(triadFeats(1:(triadTier(end)-1),:),freqs{ii})),2))/(triadTier(end)-1);
   
   topPerc(ii,1) = sum(~cellfun(@isempty,regexp(monadFeats(1:(monadTier(2)-1),:),freqs{ii})))/(monadTier(2)-1);
   topPerc(ii,2) = sum(any(~cellfun(@isempty,regexp(dyadFeats(1:(dyadTier(2)-1),:),freqs{ii})),2))/(dyadTier(2)-1);
   topPerc(ii,3) = sum(any(~cellfun(@isempty,regexp(triadFeats(1:(triadTier(2)-1),:),freqs{ii})),2))/(triadTier(2)-1);
end
%% Pre-Drinking vs. Not Drinking
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\preDrinkCombined2\','.mat');
preData = [];
notData = [];
for ii = 1:size(files,2)
   load(files{ii})
   if ~isempty(trls{1,1})
      [b,c,t] = size(psdTrls{1,1}.relPow);
      thisPow = reshape(psdTrls{1,1}.relPow,b*c,t)';
      [cmb,b,t] = size(coh{1,1}.rel);
      thisCoh = reshape(permute(coh{1,1}.rel,[2,1,3]),cmb*b,t)';
      preData = [preData;thisPow,thisCoh];
      
      [b,c,t] = size(psdTrls{1,end}.relPow);
      thisPow = reshape(psdTrls{1,end}.relPow,b*c,t)';
      [cmb,b,t] = size(coh{1,end}.rel);
      thisCoh = reshape(permute(coh{1,end}.rel,[2,1,3]),cmb*b,t)';
      notData = [notData;thisPow,thisCoh];
   end
end
%% Collate and impute
catData{1,1} = preData;
catData{1,2} = notData;
[all,each,rnd,~] = evenDataSplit(catData,17222,4306,'ADA',20);
save('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\preDrinkData.mat','all','each','rnd')
%% Build full logistics
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\preDrinkData.mat')
for n = 1:20
    disp([num2str(n),' of 20'])
    % Set up training data
    trainX = all.trainX{n};
    trainY = all.trainY{n};
    % Set up testing data
    testX = all.testX{n};
    testY = all.testY{n};
    % Build and test model on concat data
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [x{n},y{n},~,a(n)] = perfcurve(testY,prob,1);
    % Set up random testing data
    testY = rnd.allTestY{n};
    prob = predict(mdl,testX);
    [rndX{n},rndY{n},~,rndA(n)] = perfcurve(testY,prob,1);
end
%% Test drinkNot models on preDrink data
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drinkNotModel.mat','concatData')
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\preDrinkData.mat')
for ii = 1:20
    disp([num2str(ii),' of 20'])
    testX = all.testX{ii};
    testY = all.testY{ii};
    prob = predict(concatData.model{ii},testX);
    [drinkX{ii},drinkY{ii},~,drinkA(ii)] = perfcurve(testY,prob,1);
end
%%
figure
hold on
plot(mean(cat(2,x{:}),2),mean(cat(2,y{:}),2),'-k')
plot(mean(cat(2,drinkX{:}),2),mean(cat(2,drinkY{:}),2),':k')
plot(mean(cat(2,rndX{:}),2),mean(cat(2,rndY{:}),2),'--k')
plot(NaN,NaN,'color','none')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
d = distES(a,rndA);
legend({['PreDrink: ',num2str(round(mean(a),2))],['DrinkNot: ',num2str(round(mean(drinkA),2))],['Permuted: ',num2str(round(mean(rndA),2))],['d = ',num2str(round(d,2))]},'Location','se')
title('Pre-Drink vs. Other')
%% Drink times
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\mat\','.mat');
for ii = 1:length(files)
   load(files{ii},'eventTs') 
   drink{ii} = [eventTs.t{logicFind('Drink (Start)',eventTs.label,'==')},eventTs.t{logicFind('Drink (End)',eventTs.label,'==')}];
   total(ii) = eventTs.t{1,1}(end);
   norm{ii} = drink{ii}./total(ii);
end
dt = 1/max(total);
masterTime = 0:dt:1;
masterTimeline = zeros(length(files),length(masterTime));
for ii = 1:length(files)
    if ~isempty(norm{ii})
        for jj = 1:size(norm{ii},1)
            masterTimeline(ii,round(norm{ii}(jj,1)*max(total)):round(norm{ii}(jj,2)*max(total))) = 1; 
        end
    end
end
meanTimeline = mean(masterTimeline,1);
% Smooth
bins = 13;
bSize = length(masterTime)/13;
for ii = 1:bins
   smTime(ii) = mean(meanTimeline(bSize*ii-bSize+1:bSize*ii));
end
figure
plot(masterTime,smTime)
%%
smTime2 = smooth(meanTimeline,0.1,'lowess');
figure
plot(masterTime,smTime2.*100)
hold on
nanTime = masterTimeline;
nanTime(nanTime==0)=NaN;

plot(masterTime,nanTime(36,:).*0.5,'k')
box off
xlabel('Normalized Time')
ylabel('Percent of animals drinking')
%% Calculate CFC
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkNot\','.mat');
drink = [];
notDrink =[];
for ii = 1:size(files,2)
   load(files{ii}) 
   if ~isempty(trls{1,1})
       m = mean(reshape(trls{1,1}.trial,4,2000*size(trls{1,1}.trial,3)),2);
       m = repmat(m,1,2000,size(trls{1,1}.trial,3));
       s = std(reshape(trls{1,1}.trial,4,2000*size(trls{1,1}.trial,3)),[],2);
       s = repmat(s,1,2000,size(trls{1,1}.trial,3));
       drink = cat(3,drink,(trls{1,1}.trial-m)./s);
   end
   m = mean(reshape(trls{1,2}.trial,4,2000*size(trls{1,2}.trial,3)),2);
   m = repmat(m,1,2000,size(trls{1,2}.trial,3));
   s = std(reshape(trls{1,2}.trial,4,2000*size(trls{1,2}.trial,3)),[],2);
   s = repmat(s,1,2000,size(trls{1,2}.trial,3));
   notDrink = cat(3,notDrink,(trls{1,2}.trial-m)./s);
end
%%
data = reshape(drink,4,size(drink,2)*size(drink,3));
[MI,fSpace,Hvals,probVects,statsData] = gmwMI(data,400,1,100,18,size(drink,3),2000,3,5,4);
%% Try splitting and averaging
data1 = reshape(drink(:,:,1:173),4,2000*173);
data2 = reshape(drink(:,:,174:end),4,2000*173);
[MI1,fSpace1,Hvals1,probVects1,statsData1] = gmwMI(data1,400,1,100,18,173,2000,3,5,4);
[MI2,fSpace2,Hvals2,probVects2,statsData2] = gmwMI(data2,400,1,100,18,173,2000,3,5,4);
%% Try subsampling
for ii = 1:100
    subData = reshape(drink(1,:,randperm(346,100)),1,2000*100);
    [MIsub(:,:,:,:,ii),fSpace] = gmwMI(subData,400,1,100,18,100,2000,3,5,4);
end
%%
data = reshape(notDrink,4,size(notDrink,2)*size(notDrink,3));
%%
tic
[nMI,nfSpace,nHvals,nprobVects,nstatsData] = gmwMI(data(1,:),400,1,100,18,size(notDrink,3),2000,3,5,4);
toc
%%
figure
% imagesc(mean(cat(3,squeeze(MI1(:,:,1,1)'),squeeze(MI2(:,:,1,1)')),3))
% imagesc(MI(:,:,1,1)')
% imagesc(mean(MIsub(:,:,1,1,:),5)')
imagesc(mean(MIsub(:,:,1,1,:),5)'-MI(:,:,1,1)')
fSpaceH = (fSpace)*(400/(2*pi));
set(gca,'ytick',1:4:length(fSpace)), set(gca,'xtick',1:4:length(fSpace))
set(gca,'xticklabel',round(fliplr(fSpaceH([1:4:length(fSpaceH)]))',1))
set(gca,'yticklabel',round((fSpaceH([1:4:length(fSpaceH)]))',1))
set(gca,'xdir','reverse')
xlabel('Phase (Hz)')
ylabel('Amplitude (Hz)')
colormap('viridis')
box off
%% Granger Causality
load('drinkNotAllTrialsNorm.mat')
[GC,~,~,f] = condGCnPar(permute(notDrink,[2,3,1]),400);
% Time-reverse GC
[trGC,~,~,f] = condGCnPar(permute(fliplr(notDrink),[2,3,1]),400);
%%
load('')
%% Compare and Plot TRGC
ch1 = 1;
ch2 = 3;
diffGC = GC(1:100,ch1,ch2)-GC(1:100,ch2,ch1);
difftrGC = trGC(1:100,ch2,ch1)-trGC(1:100,ch1,ch2);
figure
scatter(diffGC,difftrGC,'.')
%% Permuted GC
master = [];
for ii = 1:20
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drinkGC\permtest',num2str(ii),'.mat'])
   master = cat(3,master,mGC);
end
cmbs = [nchoosek(1:4,2);fliplr(nchoosek(1:4,2))];
for ii = 1:size(cmbs,1)
    s = sort(mGC(cmbs(ii,1),cmbs(ii,2),:),'descend');
    t(cmbs(ii,1),cmbs(ii,2)) = s(1000*0.01-1);
end
%% plot
figure
chan = {'lPFC','rPFC','lS','rS'};
cmbs = [nchoosek(1:4,2);fliplr(nchoosek(1:4,2))];
for ii = 1:size(cmbs,1)
    subplot(3,4,ii)
    hold on
    plot(f,GC(:,cmbs(ii,1),cmbs(ii,2)))
    plot(f,trGC(:,cmbs(ii,1),cmbs(ii,2)),'--')
    plot(f,trGC(:,cmbs(ii,2),cmbs(ii,1)))
%     plot([0 100],[t(cmbs(ii,1),cmbs(ii,2)) t(cmbs(ii,1),cmbs(ii,2))])
    title([chan{cmbs(ii,1)},'->',chan{cmbs(ii,2)}],'FontSize',13)
    xlim([0 100])
end
%%
figure
for ii = 1:size(cmbs,1)
    subplot(3,4,ii)
    diffGC = GC(1:100,cmbs(ii,1),cmbs(ii,2))-GC(1:100,cmbs(ii,2),cmbs(ii,1));
    difftrGC = trGC(1:100,cmbs(ii,2),cmbs(ii,1))-trGC(1:100,cmbs(ii,1),cmbs(ii,2));
    scatter(diffGC,difftrGC,'.')
    hold on
    plot([0 1],[0 1])
    lsline
end
%%
figure
cmbs = nchoosek(1:4,2);
for ii = 1:6
    subplot(2,3,ii)
    hold on
    plot(f,GC(:,cmbs(ii,1),cmbs(ii,2)))
    plot(f,GC(:,cmbs(ii,2),cmbs(ii,1)))
    xlim([0 100])
    xlabel('Hz')
    ylabel('GC')
    title([num2str(cmbs(ii,1)),'-',num2str(cmbs(ii,2))])
end
%% CFC 2000 vs 346
load('drinkNotAllTrials.mat')
for ii = 1:10
    disp([num2str(ii),' of 2000'])
    grab = notDrink(:,:,randperm(1000));
    data = reshape(grab,4,size(grab,2)*size(grab,3));
    [MI,fSpace,Hvals,probVects,statsData] = gmwMI(data,400,1,100,18,size(grab,3),2000,3,5,4);
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\notDrinkCFC\notDrink_2000_',num2str(ii),'.mat'],'MI','fSpace','Hvals','probVects','statsData','-v7.3')
    clear MI fSpace Hvals probVects statsData
end
for ii = 1:10
    disp([num2str(ii),' of 346'])
    grab = notDrink(:,:,randperm(346));
    data = reshape(grab,4,size(grab,2)*size(grab,3));
    [MI,fSpace,Hvals,probVects,statsData] = gmwMI(data,400,1,100,18,size(grab,3),2000,3,5,4);
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\notDrinkCFC\notDrink_346_',num2str(ii),'.mat'],'MI','fSpace','Hvals','probVects','statsData','-v7.3')
    clear MI fSpace Hvals probVects statsData
end
%%
for ii = 1:10
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\notDrinkCFC\notDrink_346_',num2str(ii),'.mat'],'MI')
   MI346(:,:,:,:,ii) = MI;
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\notDrinkCFC\notDrink_2000_',num2str(ii),'.mat'],'MI')
   MI2000(:,:,:,:,ii) = MI;
end
%%
figure
for ii = 1:10
    imagesc(MI2000(:,:,1,2,ii)')
    fSpaceH = (fSpace)*(400/(2*pi));
    set(gca,'ytick',1:4:length(fSpace)), set(gca,'xtick',1:4:length(fSpace))
    set(gca,'xticklabel',round(fliplr(fSpaceH([1:4:length(fSpaceH)]))',1))
    set(gca,'yticklabel',round((fSpaceH([1:4:length(fSpaceH)]))',1))
    set(gca,'xdir','reverse')
    xlabel('Phase (Hz)')
    ylabel('Amplitude (Hz)')
    colormap('viridis')
    box off
    title(num2str(ii))
    pause(5)
end
%% CFC different samples
iter = 100:50:2000;
master = [];
for ii = 1:size(iter,2)
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\notDrinkCFC\notDrinkCFC_',num2str(iter(ii)),'.mat'])
    master = cat(5,master,MI);
end
%%
figure
for ii = 1:size(master,5)
    imagesc(master(:,:,4,2,ii)')
    fSpaceH = (fSpace)*(400/(2*pi));
    set(gca,'ytick',1:4:length(fSpace)), set(gca,'xtick',1:4:length(fSpace))
    set(gca,'xticklabel',round(fliplr(fSpaceH([1:4:length(fSpaceH)]))'))
    set(gca,'yticklabel',round((fSpaceH([1:4:length(fSpaceH)]))'))
    set(gca,'xdir','reverse')
    xlabel('Phase (Hz)')
    ylabel('Amplitude (Hz)')
    colormap('viridis')
    box off
    title(num2str(iter(ii)))
    colorbar
    caxis([0 18e-5])
    pause(1)
end
%%
fig = figure('pos',[750 100 1280 1024]);
cmbs = [1,1;1,2;1,3;1,4;2,1;2,2;2,3;2,4;3,1;3,2;3,3;3,4;4,1;4,2;4,3;4,4];
fSpaceH = (fSpace)*(400/(2*pi));
for ii = 1:size(master,5)
    for jj = 1:size(cmbs,1)
        subplot(4,4,jj)
        imagesc(master(:,:,cmbs(jj,1),cmbs(jj,2),ii)')
        if jj == 1
            text(130,0.9,['Samples: ',num2str(iter(ii))],'FontWeight','bold')
            colorbar('position',[0.02,0.79,0.012,0.096])
            % Add generic axes
            % Y
            annotation(fig,'line',[0.06 0.06],[0.865 0.8]);
            text(113,58,'Amp (Hz)','Rotation',90)
            % X
            annotation(fig,'line',[0.06,0.09],[0.8 0.8]);
            text(109,65,'\phi (Hz)')
            
        end
        set(gca,'ytick',1:10:length(fSpace)), set(gca,'xtick',1:10:length(fSpace))
        set(gca,'xticklabel',round(fliplr(fSpaceH([1:10:length(fSpaceH)]))'))
        set(gca,'yticklabel',round((fSpaceH([1:10:length(fSpaceH)]))'))
        title([num2str(cmbs(jj,1)),'-',num2str(cmbs(jj,2))])
        set(gca,'xdir','reverse')
        colormap('viridis')
        box off
        caxis([0 18e-5])
    end
    drawnow
    frame = getframe(fig);
    im{ii} = frame2im(frame);
end
close;
%%
filename = 'sampleCFC.gif';
for idx = 1:size(master,5)
   [A,map] = rgb2ind(im{idx},256);
   if idx == 1
       imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
   else
       imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
   end
end
%%
cfg_freqanalysis = [];
cfg_freqanalysis.output = 'fourier';
cfg_freqanalysis.method = 'mtmfft';
cfg_freqanalysis.taper = 'dpss';
cfg_freqanalysis.pad = 10;
cfg_freqanalysis.foi = 1:100;
cfg_freqanalysis.tapsmofrq = 4;
cfg_freqanalysis.keeptrials = 'yes';
cfg_freqanalysis.channel = {'lPFC','rPFC','lS','rS'};
cfg_freqanalysis.channelcmb = 'all';
% Prep data for fieldtrip format
data = [];
data.label = {'lPFC','rPFC','lS','rS'};
data.fsample = 400;
for ii = 1:size(drink,3)
    data.trial{ii} = circshift(drink(:,:,ii),randi(2000));
    data.time{ii} = (1:2000)./400;
end
addpath(genpath('C:\Users\Pythia\documents\GreenLab\code\outside\fieldtrip-master\'))
% Compute power spectrum
freq = ft_freqanalysis(cfg_freqanalysis,data);
% Compute psi
cfg_psi = [];
cfg_psi.method = 'psi';
cfg_psi.bandwidth = 2;
cfg_psi.jackknife = 'yes';
psitest = ft_connectivityanalysis(cfg_psi,freq);
%%
chan = [1,2];
figure
hold on
plot(psi2.freq,squeeze(psi2.psispctrm(chan(1),chan(2),:)),'-b')
plot(psi2.freq,squeeze(psi2.psispctrm(chan(1),chan(2),:)+psi2.psispctrmsem(chan(1),chan(2),:)),':b')
plot(psi2.freq,squeeze(psi2.psispctrm(chan(1),chan(2),:)-psi2.psispctrmsem(chan(1),chan(2),:)),':b')
xlabel('Frequency (Hz)')
ylabel(sprintf('Phase slope (%c/Hz)',char(176)))

figure
plot(psi2.freq,squeeze(psi2.psispctrm(chan(1),chan(2),:))./squeeze(psi2.psispctrmsem(chan(1),chan(2),:)),'-k')
% ylim([-10 10])
xlabel('Frequency (Hz)')
ylabel('Phase slope index')
box off