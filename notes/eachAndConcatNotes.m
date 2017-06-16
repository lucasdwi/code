%% Fix 'each' auc by applying model to zscore
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed')
load('concatData.mat')
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\300_75')
for ii = 1:20
    load([num2str(ii),'.mat'])
    concatX{ii} = concatData.rocX;
    concatY{ii} = concatData.rocY;
    concatA(ii) = concatData.auc;
    for jj = 1:12
        [predY] = cvglmnetPredict(concatData.model,zscore(eachTestX{ii,jj}),'lambda_1se','response');
        [eachX{ii,jj},eachY{ii,jj},~,eachA(ii,jj)] = perfcurve(eachTestY{ii,jj},predY,1); 
    end
    eachData.rocX = eachX;
    eachData.rocY = eachY;
    eachData.auc = eachA;
    save([num2str(ii),'fix.mat'],'concatData','eachData')
    clearvars -except ii eachTestX eachTestY
end
%% Fix 'each' from runEach5050
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(n,20);
    if iter == 0
        iter = 20;
    end
    for jj = 1:12
        if jj ~= animal
            testX = eachTestX{iter,jj};
            testY = eachTestY{iter,jj};
            [predY] = cvglmnetPredict(selfData.model,zscore(testX),'lambda_1se','response');
            [eachX{jj},eachY{jj},~,eachA(jj)] = perfcurve(testY,predY,1);
        end
    end
    eachData.rocX = eachX;
    eachData.rocY = eachY;
    eachData.auc = eachA;
    save([num2str(ii),'fix.mat'],'concatData','eachData','selfData')
end
%%
figure
plot(concatData.rocX,concatData.rocY,'-k','LineWidth',2)
hold on
for ii = 1:12
    plot(eachData.rocX{ii},eachData.rocY{ii},'-r')
end
%% Load all concat; get average for all concat->concat and concat->each
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Trials_50-50')
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Train_50-50')
for ii = 1:20
   load([num2str(ii),'.mat'])
   cA(ii) = concatData.auc;
   cX{ii} = concatData.rocX;
   cY{ii} = concatData.rocY;
   % 
   eA(ii,:) = eachData.auc;
   eX(ii,:) = eachData.rocX;
   eY(ii,:) = eachData.rocY;
end
% Get average and std of 'eachA'
eAM = mean(eA,1);
eS = std(eA,[],1);
% Get average ROC for 'each'
for ii = 1:12
    eXM(ii,:) = mean(cat(2,eX{:,ii}),2);
    eYM(ii,:) = mean(cat(2,eY{:,ii}),2);
    eXS(ii,:) = std(cat(2,eX{:,ii}),[],2);
    eYS(ii,:) = std(cat(2,eY{:,ii}),[],2);
end
%% Load all concatRand files
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Trials_50-50Rand')
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Train_50-50Rand')
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Check if roc is 50 line, if so interpolate it; store curve
   if isequal(concatData.rocX,[0;1])
      concatX{ii} = (0:1/1500:1)'; 
      concatY{ii} = (0:1/1500:1)';
   else
       concatX{ii} = concatData.rocX;
       concatY{ii} = concatData.rocY;
   end
   % Check if roc is 50 line, if so interpolate it; store curve
   if isequal(eachData.rocX,[0,1])
       eachX{ii} = (0:1/125:1)';
       eachY{ii} = (0:1/125:1)';
   else
       eachX{ii} = eachData.rocX;
       eachY{ii} = eachData.rocY;
   end
   % Store auc value
    concatA(ii) = concatData.auc;
    cEachA(ii,:) = eachData.auc;
end
% Get average ROCs
randAvgX = mean(cat(2,concatX{:}),2);
randAvgY = mean(cat(2,concatY{:}),2);
randStdX = std(cat(2,concatX{:}),[],2);
randStdY = std(cat(2,concatY{:}),[],2);
randXp = randAvgX+randStdX;
randXm = randAvgX-randStdX;
randYp = randAvgY+randStdY;
randYm = randAvgY-randStdY;
randXfill = [randXp;flipud(randXm)];
randYfill = [randYp;flipud(randYm)];
%% Compare real and permuted auc distributions using Mann-Whitney U
[dC,pC] = distES(cA,concatA);
for ii = 1:12
    [dE(ii),dP(ii)] = distES(eA(:,ii),cEachA(:,ii));
end
[~,~,pAdj] = fdr_bh(dP,0.05,'dep');
%% Plot concat ROCs
cols = distinguishable_colors(13);
figure
hold on
% for ii = 1:20
%     plot(cX{ii},cY{ii},'-k')
% end
cMX = mean(cat(2,cX{:}),2);
cMY = mean(cat(2,cY{:}),2);
cSX = std(cat(2,cX{:}),[],2);
cSY = std(cat(2,cY{:}),[],2);
cEXp = cMX+cSX;
cEXm = cMX-cSX;
cEYp = cMY+cSY;
cEYm = cMY-cSY;
cXFill = [cEXp;flipud(cEXm)];
cYFill = [cEYp;flipud(cEYm)];
plot(cMX,cMY,'Color',cols(13,:),'LineWidth',2)
fill(cXFill,cYFill,[.5 .5 .5]);
for ii = 1:12
   plot(eXM(ii,:),eYM(ii,:),'Color',cols(ii,:)) 
end
xlim([0 1])
ylim([0 1])
% Add rand line
hold on
fill(randXfill,randYfill,'r')
plot(randAvgX,randAvgY,'b')
%%
figure
for ii = 1:12
    subplot(4,3,ii)
    hold on
    for jj = 1:20
        plot(eX{jj,ii},eY{jj,ii},'-k')
    end
    plot(eXM(ii,:),eYM(ii,:),'-','Color',cols(ii,:),'LineWidth',3)
    fill(cXFill,cYFill,[.5 .5 .5]);
    xlim([0 1])
    ylim([0 1])
    %plot(mean(cat(2,cX{:}),2),mean(cat(2,cY{:}),2),'Color',cols(13,:),'LineWidth',3)
end
%% Open all concatIncrease files to determine effect of increasing training
% set while test set remains the same size and same data
for ii = 1:34
   load([num2str(ii),'.mat'])
   concatAUC(ii,:) = cell2mat(concatData.auc);
   eachAUC(:,:,ii) = cat(1,eachData.auc{:});
end
%%
mConcatA = mean(concatAUC,2);
sConcatA = std(concatAUC,[],2);
mEachA = squeeze(mean(eachAUC,1));
sEachA = squeeze(std(eachAUC,[],1));
figure
shadedErrorBar(40:20:700,mConcatA,sConcatA)
hold on
figure
hold on
for ii = 1:12
    plot(40:20:700,mEachA(ii,:))
%     shadedErrorBar(40:20:700,mEachA(ii,:),sEachA(ii,:))
end
%% Open all each files
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\500Trials_50-50\')
for ii = 1:240
    load([num2str(ii),'fix.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Get each and self auc and rocs
    eachA(animal,:,iter) = eachData.auc;
    eachA(animal,animal,iter) = selfData.auc;
    eachData.rocX{animal} = selfData.rocX;
    eachX{animal}(:,:,iter) = cat(2,eachData.rocX{:});
    eachData.rocY{animal} = selfData.rocY;
    eachY{animal}(:,:,iter) = cat(2,eachData.rocY{:});
    % Get concat auc and rocs
    concatA(animal,iter) = concatData.auc;
    concatX(animal,:,iter) = concatData.rocX;
    concatY(animal,:,iter) = concatData.rocY;
end
%%
% Average across AUCs for each
mEachA = squeeze(mean(eachA,3));
mConcatA = mean(concatA,2);
figure
imagesc([mEachA,mConcatA])
colormap('viridis')
xlabel('Test Set')
ylabel('Training Set')
set(gca,'XTick',(1:13),'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','C'});
sEachA = squeeze(std(eachA,[],3));
sConcatA = std(concatA,[],2);
figure
imagesc([sEachA,sConcatA])
colormap('viridis')
xlabel('Test Set')
ylabel('Training Set')
set(gca,'XTick',(1:13),'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','C'});
%% Plot one animal, self, each other animal, and concat to it
%%
figure
for ii = 1:12
    subplot(3,4,ii)
    selfX = mean(eachX{ii}(:,ii,:),3);
    selfY = mean(eachY{ii}(:,ii,:),3);
    plot(selfX,selfY,'b','LineWidth',2)
    hold on
    for jj = 1:12
        if jj ~= ii
            otherX{jj} = mean(eachX{jj}(:,ii,:),3);
            otherY{jj} = mean(eachY{jj}(:,ii,:),3);
            plot(otherX{jj},otherY{jj},'k')
        end
    end
    mcX = mean(cat(2,eX{:,ii}),2);
    mcY = mean(cat(2,eY{:,ii}),2);
    plot(mcX,mcY,'r','LineWidth',2)
    plot(0:1,0:1,'--k')
end
%%
figure
plot(avgX(:,1),avgY(:,1))