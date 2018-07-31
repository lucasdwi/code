load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50v2.mat')
binge = all;
bingeRnd = rnd;
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drink6000.mat','all','rnd')
drink = all;
drinkRnd = rnd;
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];
for ii = 1:20
    disp(num2str(ii))
    % Binge -> Binge
    trainX = binge.trainX{ii}(:,bingeInds);
    trainY = binge.trainY{ii};
    testX = binge.testX{ii}(:,bingeInds);
    testY = binge.testY{ii};
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [bbX{ii},bbY{ii},~,bbA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/1200:1,'UseNearest',0);
    % Binge -> Drink
    testX = drink.testX{ii}(:,drinkInds);
    testY = drink.testY{ii};
    prob = predict(mdl,testX);
    [bdX{ii},bdY{ii},~,bdA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/1200:1,'UseNearest',0);
    % Drink -> Drink
    trainX = drink.trainX{ii}(:,drinkInds);
    trainY = drink.trainY{ii};
    testX = drink.testX{ii}(:,drinkInds);
    testY = drink.testY{ii};
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [ddX{ii},ddY{ii},~,ddA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/1200:1,'UseNearest',0);
    % Drink -> Binge
    testX = binge.testX{ii}(:,bingeInds);
    testY = binge.testY{ii};
    prob = predict(mdl,testX);
    [dbX{ii},dbY{ii},~,dbA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/1200:1,'UseNearest',0);
end
%%
figure
subplot(2,2,1)
hold on
for ii = 1:20
    plot(bbX{ii},bbY{ii})
end
text(0.5,0.25,num2str(round(mean(bbA),2)))
title('Binge -> Binge')
subplot(2,2,2)
hold on
for ii = 1:20
    plot(bdX{ii},bdY{ii})
end
text(0.5,0.25,num2str(round(mean(bdA),2)))
title('Binge -> Drink')
subplot(2,2,3)
hold on
for ii = 1:20
    plot(dbX{ii},dbY{ii})
end
text(0.5,0.25,num2str(round(mean(dbA),2)))
title('Drink -> Binge')
subplot(2,2,4)
hold on
for ii = 1:20
    plot(ddX{ii},ddY{ii})
end
text(0.5,0.25,num2str(round(mean(ddA),2)))
title('Drink -> Drink')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50v2.mat')
binge = all;
bingeRnd = rnd;
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drink6000.mat','all','rnd')
drink = all;
drinkRnd = rnd;
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];
for ii = 1:20
    % Train and test gen model
    trainX = [binge.trainX{ii}(:,bingeInds);drink.trainX{ii}(:,drinkInds)];
    trainY = [binge.trainY{ii};drink.trainY{ii}];
    testX = [binge.testX{ii}(:,bingeInds);drink.testX{ii}(:,drinkInds)];
    testY = [binge.testY{ii};drink.testY{ii}];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [ggX{ii},ggY{ii},~,ggA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/2400:1,'UseNearest',0);
    % Apply gen model to binge
    testX = binge.testX{ii}(:,bingeInds);
    testY = binge.testY{ii};
    prob = predict(mdl,testX);
    [gbX{ii},gbY{ii},~,gbA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/2400:1,'UseNearest',0);
    % Apply gen model to drink
    testX = drink.testX{ii}(:,drinkInds);
    testY = drink.testY{ii};
    prob = predict(mdl,testX);
    [gdX{ii},gdY{ii},~,gdA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/2400:1,'UseNearest',0);
    % Permuted
    trainX = [binge.trainX{ii}(:,bingeInds);drink.trainX{ii}(:,drinkInds)];
    testX = [binge.testX{ii}(:,bingeInds);drink.testX{ii}(:,drinkInds)];
    trainY = [bingeRnd.allTrainY{ii};drink.trainY{ii}(randperm(size(drink.trainY{1},1)),:)];
    testY = [bingeRnd.allTestY{ii};drink.testY{ii}(randperm(size(drink.testY{1},1)),:)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [rndX{ii},rndY{ii},~,rndA(ii)] = perfcurve(testY,prob,1,'TVals',0:1/2400:1,'UseNearest',0);
end
%% Plot the six ROC curves in grid
figure
subplot(3,2,1)
plot(mean(cat(2,ddX{:}),2),mean(cat(2,ddY{:}),2),'-k')
box off
title('Drink>Drink')

subplot(3,2,2)
plot(mean(cat(2,dbX{:}),2),mean(cat(2,dbY{:}),2),'-k')
box off
title('Drink>Binge')

subplot(3,2,3)
plot(mean(cat(2,bbX{:}),2),mean(cat(2,bbY{:}),2),'-k')
box off
title('Binge>Binge')

subplot(3,2,4)
plot(mean(cat(2,bdX{:}),2),mean(cat(2,bdY{:}),2),'-k')
box off
title('Binge>Drink')

subplot(3,2,5)
plot(mean(cat(2,gdX{:}),2),mean(cat(2,gdY{:}),2),'-k')
box off
title('Gen>Drink')

subplot(3,2,6)
plot(mean(cat(2,gbX{:}),2),mean(cat(2,gbY{:}),2),'-k')
box off
title('Gen>Binge')
%% Plot gen>gen curve
figure
plot(mean(cat(2,ggX{:}),2),mean(cat(2,ggY{:}),2),'k')
hold on
plot(mean(cat(2,rndX{:}),2),mean(cat(2,rndY{:}),2),'--k')
box off
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1]);
title('Gen>Gen')
%%
figure
subplot(1,3,1)
hold on
for ii = 1:20
    plot(gbX{ii},gbY{ii})
end
title('Gen -> Binge')
subplot(1,3,2)
hold on
for ii = 1:20
   plot(ggX{ii},ggY{ii})
   plot(rndX{ii},rndY{ii},'--k')
end
title('General Consumption Model')
subplot(1,3,3)
hold on
for ii = 1:20
    plot(gdX{ii},gdY{ii})
end
title('Gen -> Drink')
%%
figure
subplot(3,3,1)
title('Binge -> Binge')
hold on
for ii = 1:20
    plot(bbX{ii},bbY{ii})
end
plot([0 1],[0 1],'--k','LineWidth',1)
text(0.8,0.1,num2str(round(mean(bbA),2)))
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)

subplot(3,3,3)
title('Binge -> Drink')
hold on
for ii = 1:20
    plot(bdX{ii},bdY{ii})
end
plot([0 1],[0 1],'--k','LineWidth',1)
text(0.8,0.1,num2str(round(mean(bdA),2)))
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)

subplot(3,3,4)
title('Gen -> Binge')
hold on
for ii = 1:20
    plot(gbX{ii},gbY{ii})
end
plot([0 1],[0 1],'--k','LineWidth',1)
text(0.8,0.1,num2str(round(mean(gbA),2)))
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)

subplot(3,3,5)
title('Gen -> Gen')
hold on
for ii = 1:20
   plot(ggX{ii},ggY{ii})
%    plot(rndX{ii},rndY{ii},'--k')
end
plot([0 1],[0 1],'--k','LineWidth',1)
text(0.8,0.1,num2str(round(mean(ggA),2)))
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)

subplot(3,3,6)
title('Gen -> Drink')
hold on
for ii = 1:20
    plot(gdX{ii},gdY{ii})
end
plot([0 1],[0 1],'--k','LineWidth',1)
text(0.8,0.1,num2str(round(mean(gdA),2)))
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)

subplot(3,3,7)
title('Drink -> Binge')
hold on
for ii = 1:20
    plot(dbX{ii},dbY{ii})
end
plot([0 1],[0 1],'--k','LineWidth',1)
text(0.8,0.1,num2str(round(mean(dbA),2)))
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)

subplot(3,3,9)
title('Drink -> Drink')
hold on
for ii = 1:20
    plot(ddX{ii},ddY{ii})
end
plot([0 1],[0 1],'--k','LineWidth',1)
text(0.8,0.1,num2str(round(mean(ddA),2)))
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
%%
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50v2.mat')
binge = all;
bingeRnd = rnd;
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drink6000.mat','all','rnd')
drink = all;
drinkRnd = rnd;
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];

for k = 1%:3
    cmbs = nchoosek(1:18,k);
    for ii = 1:20
        disp([num2str(k),': ',num2str(ii)])
        for jj = 1:size(cmbs,1)
            % Train and test gen model
            trainX = [binge.trainX{ii}(:,bingeInds);drink.trainX{ii}(:,drinkInds)];
            trainY = [binge.trainY{ii};drink.trainY{ii}];
%             testX = [binge.testX{ii}(:,bingeInds);drink.testX{ii}(:,drinkInds)];
            testXbinge = binge.testX{ii}(:,bingeInds);
            testXdrink = drink.testX{ii}(:,drinkInds);
%             testY = [binge.testY{ii};drink.testY{ii}];
            testYbinge = [binge.testY{ii}];
            testYdrink = [drink.testY{ii}];
            mdl = fitglm(trainX(:,cmbs(jj,:)),trainY,'distribution','binomial');
%             prob = predict(mdl,testX(:,cmbs(jj,:)));
%             [~,~,~,a{k}(ii,jj)] = perfcurve(testY,prob,1);
            bingeProb = predict(mdl,testXbinge(:,cmbs(jj,:)));
            [~,~,~,bingeA{k}(ii,jj)] = perfcurve(testYbinge,bingeProb,1);
            drinkProb = predict(mdl,testXdrink(:,cmbs(jj,:)));
            [~,~,~,drinkA{k}(ii,jj)] = perfcurve(testYdrink,drinkProb,1);            
        end
    end
end
%% Sort vars
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
nameVect = nameVect(bingeInds);
for k = 1%:3
   [aSort{k},sortInd{k}] = sort(mean(a{k},1)','descend');
   cmbs = nchoosek(1:18,k);
   for ii = 1:k
       vars{k}(:,ii) = nameVect(cmbs(sortInd{k},ii));
   end
end
%% Tiers
monadTier = tier(a{1}(:,sortInd{1}));
dyadTier = tier(a{2}(:,sortInd{2}));
triadTier = tier(a{3}(:,sortInd{3}));
monadFeats = vars{1};
dyadFeats = vars{2};
triadFeats = vars{3};
freqs = {'d','t','a','b','lg','hg'};
for ii = 1:size(freqs,2)
   perc(ii,1) = sum(~cellfun(@isempty,regexp(monadFeats(1:(monadTier(end)-1),:),freqs{ii})))/(monadTier(end)-1);
   perc(ii,2) = sum(any(~cellfun(@isempty,regexp(dyadFeats(1:(dyadTier(end)-1),:),freqs{ii})),2))/(dyadTier(end)-1);
   perc(ii,3) = sum(any(~cellfun(@isempty,regexp(triadFeats(1:(triadTier(end)-1),:),freqs{ii})),2))/(triadTier(end)-1);
   
   topPerc(ii,1) = sum(~cellfun(@isempty,regexp(monadFeats(1:(monadTier(2)-1),:),freqs{ii})))/(monadTier(2)-1);
   topPerc(ii,2) = sum(any(~cellfun(@isempty,regexp(dyadFeats(1:(dyadTier(2)-1),:),freqs{ii})),2))/(dyadTier(2)-1);
   topPerc(ii,3) = sum(any(~cellfun(@isempty,regexp(triadFeats(1:(triadTier(2)-1),:),freqs{ii})),2))/(triadTier(2)-1);
end
%%
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drinkNotRaw.mat')
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\bingeNotData2.mat','rawData')
bingeData{1,1} = cat(1,rawData{1,1}{:,1});
bingeData{1,2} = cat(1,rawData{1,1}{:,2});
%%
nameVect = names({'SL','SR'},{'d','t','a','b','lg','hg'});
ind = logicFind('slt',nameVect,'==');
figure
barwitherr([std(catData{1}(:,ind)),std(catData{2}(:,ind)),std(bingeData{1}(:,ind)),std(bingeData{2}(:,ind))],[mean(catData{1}(:,ind)),mean(catData{2}(:,ind)),mean(bingeData{1}(:,ind)),mean(bingeData{2}(:,ind))])
set(gca,'XTickLabel',{'Drink','Not Drink','Binge','Not Binge'})
box off
[~,p(1),~,stats{1}] = ttest2(catData{1}(:,ind),catData{2}(:,ind));
[~,p(2),~,stats{2}] = ttest2(bingeData{1}(:,ind),bingeData{2}(:,ind));
title([nameVect{ind},': ',num2str(round(stats{1}.tstat,2)),' vs. ',num2str(round(stats{2}.tstat,2))])