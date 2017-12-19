%%
inds = [1:3,5:7,10,11];
load('base.mat')
indTrainX = eachTrainX(:,inds);
indTrainY = eachTrainY(:,inds);
indTestX = eachTestX(:,inds);
indTestY = eachTestY(:,inds);
for ii = 1:20
%     genTrainX{ii} = cat(1,eachTrainX{ii,inds});
%     genTrainY{ii} = cat(1,eachTrainY{ii,inds});
    genTestX{ii} = cat(1,eachTestX{ii,inds});
    genTestY{ii} = cat(1,eachTestY{ii,inds});
end

permInd = randi(size(eachTrainX{1,1},1),16,20);
genTrainX = cell(1,20);
genTrainY = cell(1,20);
for ii = 1:20
    for jj = 1:8
       genTrainX{ii} = cat(1,genTrainX{ii},eachTrainX{ii,inds(jj)}(permInd(:,ii),:));
       genTrainY{ii} = cat(1,genTrainY{ii},eachTrainY{ii,inds(jj)}(permInd(:,ii),:));
    end
end

load('dep24.mat')
inds = [1:3,5:7,10,11];
permInd = randi(size(eachTrainX{1,1},1),16,20);
for ii = 1:20
    for jj = 1:8
        indTrainX{ii,jj} = cat(1,indTrainX{ii,jj},eachTrainX{ii,jj});
        indTrainY{ii,jj} = cat(1,indTrainY{ii,jj},eachTrainY{ii,jj});
        indTestX{ii,jj} = cat(1,indTestX{ii,jj},eachTestX{ii,jj});
        indTestY{ii,jj} = cat(1,indTestY{ii,jj},eachTestY{ii,jj});
        genTrainX{ii} = cat(1,genTrainX{ii},eachTrainX{ii,inds(jj)}(permInd(:,ii),:));
        genTrainY{ii} = cat(1,genTrainY{ii},eachTrainY{ii,inds(jj)}(permInd(:,ii),:));
    end
%     genTrainX{ii} = cat(1,genTrainX{ii},allTrainX{ii});
%     genTrainY{ii} = cat(1,genTrainY{ii},allTrainY{ii});
    genTestX{ii} = cat(1,genTestX{ii},allTestX{ii});
    genTestY{ii} = cat(1,genTestY{ii},allTestY{ii});
end
load('dep48.mat')
inds = [1:6,8,9];
permInd = randi(size(eachTrainX{1,1},1),16,20);
for ii = 1:20
    for jj = 1:8
        indTrainX{ii,jj} = cat(1,indTrainX{ii,jj},eachTrainX{ii,jj});
        indTrainY{ii,jj} = cat(1,indTrainY{ii,jj},eachTrainY{ii,jj});
        indTestX{ii,jj} = cat(1,indTestX{ii,jj},eachTestX{ii,jj});
        indTestY{ii,jj} = cat(1,indTestY{ii,jj},eachTestY{ii,jj});
        genTrainX{ii} = cat(1,genTrainX{ii},eachTrainX{ii,inds(jj)}(permInd(:,ii),:));
        genTrainY{ii} = cat(1,genTrainY{ii},eachTrainY{ii,inds(jj)}(permInd(:,ii),:));
    end
%     genTrainX{ii} = cat(1,genTrainX{ii},allTrainX{ii});
%     genTrainY{ii} = cat(1,genTrainY{ii},allTrainY{ii});
    genTestX{ii} = cat(1,genTestX{ii},allTestX{ii});
    genTestY{ii} = cat(1,genTestY{ii},allTestY{ii});
end
load('chow.mat')
inds = [1:6,8,9];
permInd = randi(size(eachTrainX{1,1},1),16,20);
for ii = 1:20
    for jj = 1:8
        indTrainX{ii,jj} = cat(1,indTrainX{ii,jj},eachTrainX{ii,jj});
        indTrainY{ii,jj} = cat(1,indTrainY{ii,jj},eachTrainY{ii,jj});
        indTestX{ii,jj} = cat(1,indTestX{ii,jj},eachTestX{ii,jj});
        indTestY{ii,jj} = cat(1,indTestY{ii,jj},eachTestY{ii,jj});
        genTrainX{ii} = cat(1,genTrainX{ii},eachTrainX{ii,inds(jj)}(permInd(:,ii),:));
        genTrainY{ii} = cat(1,genTrainY{ii},eachTrainY{ii,inds(jj)}(permInd(:,ii),:));
    end
    %     genTrainX{ii} = cat(1,genTrainX{ii},allTrainX{ii});
%     genTrainY{ii} = cat(1,genTrainY{ii},allTrainY{ii});
    genTestX{ii} = cat(1,genTestX{ii},allTestX{ii});
    genTestY{ii} = cat(1,genTestY{ii},allTestY{ii});
end
% save('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\indvGen2.mat','genTestX','genTestY','genTrainX','genTrainY','indTestX','indTestY','indTrainX','indTrainY')
%%
% files{1,1} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge','base');
% files{1,2} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge','dep24');
% files{1,3} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge','dep48');
% files{1,4} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge','chow');
% %%
% for ii = 1:size(files,2)
%     for jj = 1:size(files{1,ii},2)
%         parts = strsplit(files{ii}{jj},'_');
%         load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\mat\',parts{1},'_',parts{2},'.mat'],'eventTs')
%         sInd = logicFind(1,strcmp(eventTs.label,'Approach (Start)'),'==');
%         eInd = logicFind(1,strcmp(eventTs.label,'Approach (End)'),'==');
%         app{ii,jj} = eventTs.t{1,eInd}-eventTs.t{1,sInd};
%     end
% end
% %%
% base = cat(1,app{1,:});
% dep24 = cat(1,app{2,:});
% dep48 = cat(1,app{3,:});
% chow = cat(1,app{4,:});
%% Find approach times associated with binge sessions of all pre-feeding windows
files{1,1} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','base');
files{1,2} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','dep24');
files{1,3} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','dep48');
files{1,4} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','chow');
for ii = 1:size(files,2)
    for jj = 1:size(files{1,ii},2)
        load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined\',files{ii}{jj}],'trls')
        for m = 1:61
            bStart = trls{m}.sampleinfo(:,2)./400;
            parts = strsplit(files{ii}{jj},'_');
            load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\mat\',parts{1},'_',parts{2}],'eventTs')
            sInd = logicFind(1,strcmp(eventTs.label,'Approach (Start)'),'==');
            eInd = logicFind(1,strcmp(eventTs.label,'Approach (End)'),'==');
            for k = 1:size(bStart,1)
                ind = logicFind(1,eventTs.t{1,eInd}+40>bStart(k) & eventTs.t{1,eInd}<bStart(k),'==');
                if ~isempty(ind)
                    app{ii}{m,jj}(k) = eventTs.t{1,eInd}(ind(end))-eventTs.t{1,sInd}(ind(end));
                else
                    app{ii}{m,jj}(k) = NaN;
                end
            end
        end
    end
end
%% Find the percent of trials coming from binge sessions with approaches that overlap with pre-feeding windows
eoi = 5:65;
    for ii = 1:61
        this = [];
        for jj = 1:4 
            this = [this,(cat(2,app{jj}{ii,:})<eoi(ii) & cat(2,app{jj}{ii,:})>eoi(ii)-5)];
        end
        perc(ii) = sum(this)/numel(this);
    end
figure
plot(2.5:62.5,perc.*100,'k')
set(gca,'xtick',2.5:10:62.5,'YTick',0:10:30)
xlabel('Time before Feeding (sec)')
ylabel('Percent of Data from Approach Behavior (%)')
box off
%%
% combineApp = [];
% for ii = 1:12
%     combineApp = cat(2,combineApp,unique(cat(2,app{:,ii})));
% end
% nonNaNApp = combineApp(~isnan(combineApp));
% eoi  = 5:40;
% for ii = 1:36
%     num(ii) = sum(nonNaNApp < eoi(ii) & nonNaNApp > eoi(ii)-5);
% end
% figure
% plot(2.5:37.5,num)
%% Plotting baseline binge size
% Binge - Rest
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\baseBingeSize\baseBingeSizeB-R2.mat')
doubleHist(real{1,1}.err,perm{1,1}.err(1:1000),'Main','Predicting Baseline Binge Size: Binge - Rest','xlab','Mean Absolute Error (gm)');
% d = distES(real{1,1}.err,perm{1,1}.err(1:1000));
% figure
% histogram(real{1,1}.err,'FaceColor','k','Normalization','probability','BinWidth',.1,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm{1,1}.err(1:1000),'FaceColor','w','Normalization','probability','BinWidth',.1,'FaceAlpha',1)
% xlim([1 8])
% box off
% title('Predicting Baseline Binge Size: Binge - Rest')
% xlabel('Mean Absolute Error (gm)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real{1,1}.err),2)),'\pm',num2str(round(std(real{1,1}.err),2)),' gm'],['Permuted: \mu = ',num2str(round(mean(perm{1,1}.err(1:1000)),2)),'\pm',num2str(round(std(perm{1,1}.err(1:1000)),2)),' gm']})
% text(4,.11,['d = ',num2str(d)])
% Binge
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\baseBingeSize\baseBingeSizeB2.mat')
doubleHist(real{1,1}.err,perm{1,1}.err(1:1000),'Main','Predicting Baseline Binge Size: Binge','xlab','Mean Absolute Error (gm)')
% d = distES(real{1,1}.err,perm{1,1}.err(1:1000));
% figure
% histogram(real{1,1}.err,'FaceColor','k','Normalization','probability','BinWidth',.1,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm{1,1}.err(1:1000),'FaceColor','w','Normalization','probability','BinWidth',.1,'FaceAlpha',1)
% xlim([1 8])
% box off
% title('Predicting Baseline Binge Size: Binge')
% xlabel('Mean Absolute Error (gm)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real{1,1}.err),2)),'\pm',num2str(round(std(real{1,1}.err),2)),' gm'],['Permuted: \mu = ',num2str(round(mean(perm{1,1}.err(1:1000)),2)),'\pm',num2str(round(std(perm{1,1}.err(1:1000)),2)),' gm']})
% text(4,.11,['d = ',num2str(d)])
% Rest
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\baseBingeSize\baseBingeSizeR2.mat')
doubleHist(real{1,1}.err,perm{1,1}.err(1:1000),'Main','Predicting Baseline Binge Size: Rest','xlab','Mean Absolute Error (gm)')
% d = distES(real{1,1}.err,perm{1,1}.err(1:1000));
% figure
% histogram(real{1,1}.err,'FaceColor','k','Normalization','probability','BinWidth',.1,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm{1,1}.err(1:1000),'FaceColor','w','Normalization','probability','BinWidth',.1,'FaceAlpha',1)
% xlim([1 8])
% box off
% title('Predicting Baseline Binge Size: Rest')
% xlabel('Mean Absolute Error (gm)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real{1,1}.err),2)),'\pm',num2str(round(std(real{1,1}.err),2)),' gm'],['Permuted: \mu = ',num2str(round(mean(perm{1,1}.err(1:1000)),2)),'\pm',num2str(round(std(perm{1,1}.err(1:1000)),2)),' gm']})
% text(4,.11,['d = ',num2str(d)])
%% Plotting baseline binge size change
% Binge - Rest
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\bingeSizeChange\bingeSizeChangeB-R2.mat')
doubleHist(real{1,1}.err,perm{1,1}.err(1:1000),'Main','Predicting Binge Size Change: Binge - Rest','xlab','Mean Absolute Error (%)');
% d = distES(real{1,1}.err,perm{1,1}.err(1:1000));
% figure
% histogram(real{1,1}.err,'FaceColor','k','Normalization','probability','BinWidth',.05,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm{1,1}.err(1:1000),'FaceColor','w','Normalization','probability','BinWidth',.05,'FaceAlpha',1)
% xlim([0.5 2.5])
% box off
% title('Predicting Baseline Binge Size Change: Binge - Rest')
% xlabel('Mean Absolute Error (gm)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real{1,1}.err),2)),'\pm',num2str(round(std(real{1,1}.err),2)),' gm'],['Permuted: \mu = ',num2str(round(mean(perm{1,1}.err(1:1000)),2)),'\pm',num2str(round(std(perm{1,1}.err(1:1000)),2)),' gm']})
% text(1.45,.27,['d = ',num2str(d)])
% Binge
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\bingeSizeChange\bingeSizeChangeB2.mat')
doubleHist(real{1,1}.err,perm{1,1}.err(1:1000),'Main','Predicting Binge Size Change: Binge','xlab','Mean Absolute Error (%)');
% d = distES(real{1,1}.err,perm{1,1}.err(1:1000));
% figure
% histogram(real{1,1}.err,'FaceColor','k','Normalization','probability','BinWidth',.05,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm{1,1}.err(1:1000),'FaceColor','w','Normalization','probability','BinWidth',.05,'FaceAlpha',1)
% xlim([0.5 2.5])
% box off
% title('Predicting Baseline Binge Size Change: Binge')
% xlabel('Mean Absolute Error (gm)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real{1,1}.err),2)),'\pm',num2str(round(std(real{1,1}.err),2)),' gm'],['Permuted: \mu = ',num2str(round(mean(perm{1,1}.err(1:1000)),2)),'\pm',num2str(round(std(perm{1,1}.err(1:1000)),2)),' gm']})
% text(1.45,.22,['d = ',num2str(d)])
% Rest
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\bingeSizeChange\bingeSizeChangeR2.mat')
doubleHist(real{1,1}.err,perm{1,1}.err(1:1000),'Main','Predicting Binge Size Change: Rest','xlab','Mean Absolute Error (%)');
% d = distES(real{1,1}.err,perm{1,1}.err(1:1000));
% figure
% histogram(real{1,1}.err,'FaceColor','k','Normalization','probability','BinWidth',.05,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm{1,1}.err(1:1000),'FaceColor','w','Normalization','probability','BinWidth',.05,'FaceAlpha',1)
% xlim([0.5 2.5])
% box off
% title('Predicting Baseline Binge Size Change: Rest')
% xlabel('Mean Absolute Error (gm)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real{1,1}.err),2)),'\pm',num2str(round(std(real{1,1}.err),2)),' gm'],['Permuted: \mu = ',num2str(round(mean(perm{1,1}.err(1:1000)),2)),'\pm',num2str(round(std(perm{1,1}.err(1:1000)),2)),' gm']})
% text(1.45,.21,['d = ',num2str(d)])
%% Plotting palatability
% Binge - Rest
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\palatability\palatB-R3.mat')
doubleHist(1-real.err,1-perm.err(1:1000),'Main','Predicting Palatability: Binge-Rest','xlab','Accuracy (%)')
% real.acc = (1-real.err)*100;
% perm.acc = (1-perm.err)*100;
% d = distES(real.acc,perm.acc(1:1000));
% figure
% histogram(real.acc,'FaceColor','k','Normalization','probability','BinWidth',2,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm.acc(1:1000),'FaceColor','w','Normalization','probability','BinWidth',2,'FaceAlpha',1)
% xlim([0 90])
% box off
% title('Predicting Palatability: Binge - Rest')
% xlabel('Accuracy (%)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real.acc),2)),'\pm',num2str(round(std(real.acc),2)),'%'],['Permuted: \mu = ',num2str(round(mean(perm.acc(1:1000)),2)),'\pm',num2str(round(std(perm.acc(1:1000)),2)),'%']},'Location','northwest')
% text(12,.067,['d = ',num2str(d)])
% Binge
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\palatability\palatB2.mat')
doubleHist(1-real.err,1-perm.err(1:1000),'Main','Predicting Palatability: Binge','xlab','Accuracy (%)')
% real.acc = (1-real.err)*100;
% perm.acc = (1-perm.err)*100;
% d = distES(real.acc,perm.acc(1:1000));
% figure
% histogram(real.acc,'FaceColor','k','Normalization','probability','BinWidth',2,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm.acc(1:1000),'FaceColor','w','Normalization','probability','BinWidth',2,'FaceAlpha',1)
% xlim([0 90])
% box off
% title('Predicting Palatability: Binge')
% xlabel('Accuracy (%)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real.acc),2)),'\pm',num2str(round(std(real.acc),2)),'%'],['Permuted: \mu = ',num2str(round(mean(perm.acc(1:1000)),2)),'\pm',num2str(round(std(perm.acc(1:1000)),2)),'%']},'Location','northwest')
% text(12,.056,['d = ',num2str(d)])
% Rest
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\palatability\palatR2.mat')
doubleHist(1-real.err,1-perm.err(1:1000),'Main','Predicting Palatability: Rest','xlab','Accuracy (%)')
% real.acc = (1-real.err)*100;
% perm.acc = (1-perm.err)*100;
% d = distES(real.acc,perm.acc(1:1000));
% figure
% histogram(real.acc,'FaceColor','k','Normalization','probability','BinWidth',2,'FaceAlpha',1,'EdgeColor','w')
% hold on
% histogram(perm.acc(1:1000),'FaceColor','w','Normalization','probability','BinWidth',2,'FaceAlpha',1)
% xlim([0 90])
% box off
% title('Predicting Palatability: Rest')
% xlabel('Accuracy (%)')
% ylabel('Model Frequency')
% legend({['Real: \mu = ',num2str(round(mean(real.acc),2)),'\pm',num2str(round(std(real.acc),2)),'%'],['Permuted: \mu = ',num2str(round(mean(perm.acc(1:1000)),2)),'\pm',num2str(round(std(perm.acc(1:1000)),2)),'%']},'Location','northeast')
% text(40,.076,['d = ',num2str(d)])
%% Prep binge vs. not binge trial data
[data,~] = collateData('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge\',{'base','dep24','dep48','chow'},{'pow','coh'},'trl');
for ii = 1:size(data,2)
    for k = 1:size(data{1,ii},2)
        trialDat{k,ii} = cat(1,data{1,ii}{:,k});
    end
end
%% Prep average data
[data,~] = collateData('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge\',{'base','dep24','dep48','chow'},{'pow','coh','corr'},'avg','');
for ii = 1:4
    for k = 1:2
        avgDat{k,ii} = cat(1,data{1,ii}{:,k});
    end
end
%% Find potential noise features
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\bingeNotData.mat')
%%
% Load voracity data
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\4conditionBingeSize.mat')
% Correlate voracity with all features
vorDiff = [(voracity(:,3)-voracity(:,1))./voracity(:,1);(voracity(:,2)-voracity(:,1))./voracity(:,1)];
% Get non-NaN indices
inds24 = logicFind(1,~isnan(voracity(:,2)),'==');
inds48 = logicFind(1,~isnan(voracity(:,3)),'==');
% Do same subtractions for features
dep48Feats = (avgDat{1,3}-avgDat{2,3});
dep24Feats = (avgDat{1,2}-avgDat{2,2});
baseFeats = (avgDat{1,1}-avgDat{2,1});

feats = [(dep48Feats-baseFeats(inds48,:))./baseFeats(inds48,:);(dep24Feats-baseFeats(inds24,:))./baseFeats(inds24,:)];
for fi = 1:size(feats,2)
    [thisR,thisP] = corrcoef(vorDiff(~isnan(vorDiff)),feats(:,fi));
    r(fi) = thisR(1,2)^2;
    p(fi) = thisP(1,2);
end
% Get indices of significant p-values to exclude those features
pInds = logicFind(0.05,p,'<=');
%% Plot features with significant correlations with changes in voracity
% load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\bingeNotData.mat')
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
% Add worst p-value to list
for ii = 1:length(pInds)
    figure
    % Convert features to gm/ms and %
    scatter(1000.*vorDiff(~isnan(vorDiff)),100.*feats(:,pInds(ii)),'ok','Filled')
    lsline
    thiscorr{ii} = fitlm(1000.*vorDiff(~isnan(vorDiff)),feats(:,pInds(ii)));
    title(nameVect{pInds(ii)})
    text(0,0,['R^2 = ',num2str(round(thiscorr{ii}.Rsquared.Ordinary,2)),newline,'p = ',num2str(round(thiscorr{ii}.Coefficients.pValue(2),3))])
    xlabel('Change in Voracity (gm/ms)')
    ylabel('Percent Change in Feature')
end
%% Load all concat files and process
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\concat\')
% cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50')
% cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\concatOver')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\Base1')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\57\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBinge\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    ccA(ii) = concatData.auc;
    ccX{ii} = concatData.rocX;
    ccY{ii} = concatData.rocY;
    %
    ceA(ii,:) = eachData.auc;
    ceX(ii,:) = eachData.rocX;
    ceY(ii,:) = eachData.rocY;
end
% Get average and std of 'eachA'
ceAM = mean(ceA,1);
ceS = std(ceA,[],1);
% Get average ROC for 'each'
for ii = 1:12
    ceXM(ii,:) = mean(cat(2,ceX{:,ii}),2);
    ceYM(ii,:) = mean(cat(2,ceY{:,ii}),2);
    ceXS(ii,:) = std(cat(2,ceX{:,ii}),[],2);
    ceYS(ii,:) = std(cat(2,ceY{:,ii}),[],2);
end
% Get average ROC for concat-concat models
ccMX = mean(cat(2,ccX{:}),2);
ccMY = mean(cat(2,ccY{:}),2);
% Get average 'fill' for concat-concat models
[ccXfill,ccYfill] = avgFill(cat(2,ccX{:}),cat(2,ccY{:}),2,1);
%% Load all concatRand files
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\concatRand\')
% cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50Rand')
% cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\concatOverRand')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\base1Rand\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\57Rand\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBingeRand\')
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Check if roc is 50 line, if so interpolate it; store curve
   if isequal(concatData.rocX,[0;1])
      ccRandX{ii} = (0:1/1200:1)'; 
      ccRandY{ii} = (0:1/1200:1)';
   else
       ccRandX{ii} = concatData.rocX;
       ccRandY{ii} = concatData.rocY;
   end
   % Check if roc is 50 line, if so interpolate it; store curve
   for jj = 1:12
       if isequal(eachData.rocX{jj},[0;1])
           ceRandX{ii}{jj} = (0:1/125:1)';
           ceRandY{ii}{jj} = (0:1/125:1)';
       else
           ceRandX{ii}{jj} = eachData.rocX;
           ceRandY{ii}{jj} = eachData.rocY;
       end
   end
   % Store auc value
    ccRandA(ii) = concatData.auc;
    ceRandA(ii,:) = eachData.auc;
end
% Get average ROC for concat-concat random models
ccRandMX = mean(cat(2,ccRandX{:}),2);
ccRandMY = mean(cat(2,ccRandY{:}),2);
% Get average 'fill' for concat-concat random models
[ccRandXfill,ccRandYfill] = avgFill(cat(2,ccRandX{:}),cat(2,ccRandY{:}),2,1);
%%
d = distES(ccA,ccRandA);
figure
hold on
fill(ccXfill,ccYfill,[.8 .8 .8])
h(1) = plot(ccMX,ccMY,'-k');
fill(ccRandXfill,ccRandYfill,[.8 .8 .8])
h(2) = plot(ccRandMX,ccRandMY,'-k');
h(3) = plot(NaN,NaN,'Color','none');
legend(h,{['Real: ',num2str(round(mean(ccA),2))],['Permuted: ',num2str(round(mean(ccRandA),2))],['d = ',num2str(round(d,2))]},'Location','southeast')
ylim([0 1])
xlim([0 1])
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Normal Oversample')
%% Load all 'each' files
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\each\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\base1\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\57\')
for ii = 1:240
    load([num2str(ii),'fix.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Get each and self auc and rocs
    eachData.auc(animal) = selfData.auc;
    eeA(animal,:,iter) = eachData.auc;
    eachData.rocX{animal} = selfData.rocX;
    eeX{animal}(:,:,iter) = cat(2,eachData.rocX{:});
    eachData.rocY{animal} = selfData.rocY;
    eeY{animal}(:,:,iter) = cat(2,eachData.rocY{:});
    % Get concat auc and rocs
    ecA(animal,iter) = concatData.auc;
    ecX(animal,:,iter) = concatData.rocX;
    ecY(animal,:,iter) = concatData.rocY;
end
%% Load all 'each' random files 
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\eachRand\\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\57Rand\')
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Add self data to each
    eachData.rocX{animal} = selfData.rocX;
    eachData.rocY{animal} = selfData.rocY;
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(concatData.rocX,[0;1])
        ecRandX{animal,iter} = (0:1/1500:1)';
        ecRandY{animal,iter} = (0:1/1500:1)';
    else
        ecRandX{animal,iter} = concatData.rocX;
        ecRandY{animal,iter} = concatData.rocY;
    end
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(eachData.rocX{animal},[0,1])
        eeRandX{animal,iter} = (0:1/125:1)';
        eeRandY{animal,iter} = (0:1/125:1)';
    else
        eeRandX{animal,iter} = eachData.rocX;
        eeRandY{animal,iter} = eachData.rocY;
    end
    % Get each and self AUC
    eachData.auc(animal) = selfData.auc;
    eeRandA(animal,:,iter) = eachData.auc;
    % Get concat AUC
    ecRandA(ii) = concatData.auc;
end
%% Plot 'selfs' - diag
cols = distinguishable_colors(12);
figure
hold on
for ii = 1:12
   selfXm(ii,:) = mean(eeX{1,ii}(:,ii,:),3);
   selfYm(ii,:) = mean(eeY{1,ii}(:,ii,:),3);
   plot(mean(eeX{1,ii}(:,ii,:),3),mean(eeY{1,ii}(:,ii,:),3),'Color',cols(ii,:)) 
end
plot(mean(selfXm,1),mean(selfYm,1),'-k','LineWidth',2)
plot([0 1],[0 1],'--k')
title('Each to Self: Diagonal')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
text(0.8,0.2,['AUC: ',num2str(round(mean(diag(mean(eeA,3))),2))])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
%% Plot inds to others - column
cols = distinguishable_colors(12);
mX = squeeze(mean(eeX{1,1},3));
mY = squeeze(mean(eeY{1,1},3));
% Get average curve
otherXM = mean(mX,2);
otherYM = mean(mY,2);
figure
hold on
for ii = 1:12
    plot(mX(:,ii),mY(:,ii),'Color',cols(ii,:))
end
% Plot average
plot(otherXM,otherYM,'-k','LineWidth',2)
plot([0 1],[0 1],'--k')
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
text(0.8,0.2,['AUC: ',num2str(round(mean(mean(eeA(:,1,:),3),1),2))])
title('Each to Each: Column')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
%% Plot all ESs
% Convert AUCs to ESs
for ii = 1:12
    for jj = 1:12
        [eeES(ii,jj),~] = distES(squeeze(eeA(ii,jj,:)),squeeze(eeRandA(ii,jj,:))); 
    end
end
figure
pcolor(padarray(eeES,[1,1],'post'))
set(gca,'XTick',1.5:12.5,'XTickLabel',1:12,'YTick',1.5:12.5,'YTickLabel',1:12)
title('Effect Size')
xlabel('Training Set')
ylabel('Test Set')
colormap('viridis')
cb = colorbar;
ylabel(cb,'Cohen''s d')
set(cb,'FontName','Arial','FontSize',12)
%% Plot all AUCs
% Average
figure
pcolor(padarray(mean(eeA,3),[1,1],'post'))
set(gca,'XTick',1.5:12.5,'XTickLabel',1:12,'YTick',1.5:12.5,'YTickLabel',1:12)
colormap('viridis')
cb = colorbar;
ylabel(cb,'AUC')
set(cb,'FontName','Arial','FontSize',12)
title('AUC: Average')
ylabel('Training Set')
xlabel('Test Set')
% Standard deviation
figure
pcolor(padarray(std(eeA,[],3),[1,1],'post'))
set(gca,'XTick',1.5:12.5,'XTickLabel',1:12,'YTick',1.5:12.5,'YTickLabel',1:12)
colormap('viridis')
cb = colorbar;
ylabel(cb,'AUC')
set(cb,'FontName','Arial','FontSize',12)
title('AUC: Standard Deviation')
ylabel('Training Set')
xlabel('Test Set')
%% Plot concat-concat; concat-each
% Get effect size of concat-concat
[d,~] = distES(ccA,ccRandA);
figure
hold on
% Concat-each
for ii = 1:12
   plot(ceXM(ii,:),ceYM(ii,:),'Color',cols(ii,:)) 
end
% Concat-concat fill
fill(ccXfill,ccYfill,[0.8 0.8 0.8])
% Concat-concat average
plot(ccMX,ccMY,'-k','LineWidth',2)
% Concat-concat fill permuted
fill(ccRandXfill,ccRandYfill,[0.8 0.8 0.8])
% Concat-concat average permuted
plot(ccRandMX,ccRandMY,'--k','LineWidth',2)
xlim([0 1])
ylim([0 1])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
title('Concat vs. Concat-Each-Perm')
xlabel('False Positive Rate')
ylabel('True Positive Rate')

text(0.7,0.2,['d = ',num2str(round(d,2))])
text(0.7,0.15,['AUC = ',num2str(round(mean(ccA),2))])
%% Plot example of each-self vs. concat-each
% Get average each-each
for ii = 1:12
   eeXM(:,:,ii) = mean(eeX{ii},3);
   eeYM(:,:,ii) = mean(eeY{ii},3);
end
% Get AUC diffs
diffAUC = diag(mean(eeA,3))-mean(ceA,1)';
% Plot both ROCs and area
figure
hold on
fill([eeXM(:,1,1);flipud(ceXM(1,:)')],[eeYM(:,1,1);flipud(ceYM(1,:)')],[0.8 0.8 0.8]);
plot(eeXM(:,1,1),eeYM(:,1,1),'-k','LineWidth',2)
plot(ceXM(1,:),ceYM(1,:),'--k','LineWidth',2)
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
title('Concat vs. Self: Example')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend({['AUC Diff = ',num2str(round(diffAUC(1),2))],'Self','Concat'})
% Plot distribution of auc differences
figure
scatter(ones(1,12),diffAUC,400,'.k')
set(gca,'XTick',[],'XTickLabel',{})
title('AUC Differences: Self-Concat')
ylabel('AUC Difference')
%% Plot SVM data
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\svmNew.mat')
ccSVMa = concatA;
ccSVMx = mean(cat(2,concatX{:}),2);
ccSVMy = mean(cat(2,concatY{:}),2);
[ccSVMxFill,ccSVMyFill] = avgFill(cat(2,concatX{:}),cat(2,concatY{:}),2);
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\svmRandNew.mat')
ccSVMaRand = concatA;
ccSVMxRand = mean(cat(2,concatX{:}),2);
ccSVMyRand = mean(cat(2,concatY{:}),2);
[ccSVMxFillRand,ccSVMyFillRand] = avgFill(cat(2,concatX{:}),cat(2,concatY{:}),2);
% Get effect size
[svmES,~] = distES(ccSVMa,ccSVMaRand);
%%
figure
hold on
% Add concat-concat Lasso
% fill(ccXfill,ccYfill,[0.8 0.8 0.8])
% h1 = plot(ccMX,ccMY,':k','LineWidth',2);
% Then SVM
fill(ccSVMxFill,ccSVMyFill,[0.8 0.8 0.8])
h2 = plot(ccSVMx,ccSVMy,'-k','LineWidth',2);
fill(ccSVMxFillRand,ccSVMyFillRand,[0.8 0.8 0.8])
h3 = plot(ccSVMxRand,ccSVMyRand,'--k','LineWidth',2);
legend([h2,h3],{'SVM','Permuted'},'Location','southeast')
xlim([0 1])
ylim([0 1])
text(0.8,0.3,['d = ',num2str(round(svmES,2))])
text(0.8,0.25,['AUC = ',num2str(round(mean(ccSVMa),2))])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
xlabel('False Postive Rate')
ylabel('True Positive Rate')
title('Concat: SVM')
%% Compare real and permuted auc distributions using Mann-Whitney U
[dC,pC] = distES(ccA,ccRandA);
for ii = 1:12
    [dE(ii),pE(ii)] = distES(ceA(:,ii),ceRandA(:,ii));
end
[~,~,pAdj] = fdr_bh([pC,pE],0.05,'dep');
%% Build and test straight logistic models
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\concat50-50_500TrainNew.mat')
clear concatX concatY concatA
for ii = 1:20
    mdl{ii} = fitglm(allTrainX{ii},allTrainY{ii},'distribution','binomial');
    prob(:,ii) = predict(mdl{ii},allTestX{ii});
    [concatX(:,ii),concatY(:,ii),~,concatA(ii)] = perfcurve(allTestY{ii},prob(:,ii),1);
end
ccLogA = concatA;
ccLogMX = mean(concatX,2);
ccLogMY = mean(concatY,2);
[ccLogXfill,ccLogYfill] = avgFill(concatX,concatY,2);
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\concat50-50_500TrainNewRand.mat')
for ii = 1:20
    mdlRand{ii} = fitglm(allTrainX{ii},allTrainY{ii},'distribution','binomial');
    prob(:,ii) = predict(mdlRand{ii},allTestX{ii});
    [concatRandX(:,ii),concatRandY(:,ii),~,concatRandA(ii)] = perfcurve(allTestY{ii},prob(:,ii),1);
end
ccLogRandA = concatRandA;
ccLogMXrand = mean(concatRandX,2);
ccLogMYrand = mean(concatRandY,2);
[ccLogXrandFill,ccLogYrandFill] = avgFill(concatRandX,concatRandY,2);
[ccES,~] = distES(ccLogA,ccLogRandA);
%%  Plot concat-concat model with permuted data
figure
hold on
% Add concat-concat Lasso
% fill(ccXfill,ccYfill,[0.8 0.8 0.8])
% h1 = plot(ccMX,ccMY,':k','LineWidth',2);
% Then full logisitic
fill(ccLogXfill,ccLogYfill,[0.8 0.8 0.8])
h2 = plot(ccLogMX,ccLogMY,'-k','LineWidth',2);
fill(ccLogXrandFill,ccLogYrandFill,[0.8 0.8 0.8])
h3 = plot(ccLogMXrand,ccLogMYrand,'--k','LineWidth',2);
legend([h2,h3],{'Logistic','Permuted'},'Location','southeast')
xlim([0 1])
ylim([0 1])
xlabel('False Positive Rate')
ylabel('True Positive Rate')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Concat-Concat: Logistic')
text(0.8,.3,['d = ',num2str(round(ccES,2))])
text(0.8,.25,['AUC = ',num2str(round(mean(concatA),2))])
%% Test models made from concat on other conditions
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\bingeNotData.mat','pInds')
% Get usable indices
inds = 1:60;
inds = inds(~ismember(inds,pInds));
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\dep24TestNew.mat')
for ii = 1:20
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50\',num2str(ii),'.mat'])
   testX = [];
   testY = [];
   testX = allTestX{ii}(:,inds);
   testY = allTestY{ii};
   prob = cvglmnetPredict(concatData.model,zscore(testX),'lambda_1se','response');
   [dep24X{ii},dep24Y{ii},~,dep24A(ii)] = perfcurve(testY,prob,1);
   dep24ConcatES = distES(ccA,dep24A);
end
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\dep48TestNew.mat')
for ii = 1:20
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50\',num2str(ii),'.mat'])
   testX = [];
   testY = [];
   testX = allTestX{ii}(:,inds);
   testY = allTestY{ii};
   prob = cvglmnetPredict(concatData.model,zscore(testX),'lambda_1se','response');
   [dep48X{ii},dep48Y{ii},~,dep48A(ii)] = perfcurve(testY,prob,1);
   dep48ConcatES = distES(ccA,dep48A);
end
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\chowTestNew.mat')
for ii = 1:20
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50\',num2str(ii),'.mat'])
   testX = [];
   testY = [];
   testX = allTestX{ii}(:,inds);
   testY = allTestY{ii};
   prob = cvglmnetPredict(concatData.model,zscore(testX),'lambda_1se','response');
   [chowX{ii},chowY{ii},~,chowA(ii)] = perfcurve(testY,prob,1);
   chowConcatES = distES(ccA,chowA);
end
%% Prepare for plotting
dep24XM = mean(cat(2,dep24X{:}),2);
dep24YM = mean(cat(2,dep24Y{:}),2);
[dep24Xfill,dep24Yfill] = avgFill(cat(2,dep24X{:}),cat(2,dep24Y{:}),2);
dep48XM = mean(cat(2,dep48X{:}),2);
dep48YM = mean(cat(2,dep48Y{:}),2);
[dep48Xfill,dep48Yfill] = avgFill(cat(2,dep48X{:}),cat(2,dep48Y{:}),2);
chowXM = mean(cat(2,chowX{:}),2);
chowYM = mean(cat(2,chowY{:}),2);
[chowXfill,chowYfill] = avgFill(cat(2,chowX{:}),cat(2,chowY{:}),2);
%% Plot
figure
hold on
% Add concat-concat line
plot(ccMX,ccMY,'-k','LineWidth',2)
% fill(dep24Xfill,dep24Yfill,'b','FaceAlpha',0.5)
plot(dep24XM,dep24YM,'-k','LineWidth',1)
% fill(dep48Xfill,dep48Yfill,'r','FaceAlpha',0.5)
plot(dep48XM,dep48YM,'--k','LineWidth',1.5)
% fill(chowXfill,chowYfill,'y','FaceAlpha',0.5)
plot(chowXM,chowYM,':k','LineWidth',1.5)
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend({['Base: ',num2str(round(mean(ccA),2))],['Dep24: ',num2str(round(mean(dep24A),2))],['Dep48: ',num2str(round(mean(dep48A),2))],['Chow: ',num2str(round(mean(chowA),2))]},'Location','southeast')
title('Baseline Model across Conditions')
ylim([0 1])
xlim([0 1])
%% Test models made from each on other conditions
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\dep24TestNew.mat')
dep24TestX = eachTestX;
dep24TestY = eachTestY;
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\dep48TestNew.mat')
dep48TestX = eachTestX;
dep48TestY = eachTestY;
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\chowTestNew.mat')
chowTestX = eachTestX;
chowTestY = eachTestY;
inds = 1:60;
inds = inds(~ismember(inds,pInds));
for ii = 1:240
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\each\',num2str(ii),'.mat'])
   animal = ceil(ii/20);
   iter = rem(ii,20);
   if iter == 0
       iter = 20;
   end
   testX = [];
   testY = [];
   testX = dep24TestX{iter,animal}(:,inds);
   testY = dep24TestY{iter,animal};
   prob = cvglmnetPredict(selfData.model,zscore(testX),'lambda_1se','response');
   [dep24XEach{iter,animal},dep24YEach{iter,animal},~,dep24AEach(iter,animal)] = perfcurve(testY,prob,1);
   if animal <= 9
       testX = [];
       testY = [];
       testX = dep48TestX{iter,animal}(:,inds);
       testY = dep48TestY{iter,animal};
       prob = cvglmnetPredict(selfData.model,zscore(testX),'lambda_1se','response');
       [dep48XEach{iter,animal},dep48YEach{iter,animal},~,dep48AEach(iter,animal)] = perfcurve(testY,prob,1);
       testX = [];
       testY = [];
       testX = chowTestX{iter,animal}(:,inds);
       testY = chowTestY{iter,animal};
       prob = cvglmnetPredict(selfData.model,zscore(testX),'lambda_1se','response');
       [chowXEach{iter,animal},chowYEach{iter,animal},~,chowAEach(iter,animal)] = perfcurve(testY,prob,1);
   end
end
%% Prepare for plotting
[dep24xFillEach,dep24yFillEach] = avgFill(cat(2,dep24XEach{:}),cat(2,dep24YEach{:}),2);
dep24XMEach = mean(cat(2,dep24XEach{:}),2);
dep24YMEach = mean(cat(2,dep24YEach{:}),2);
[dep24ES,~] = distES(dep24A,reshape(dep24AEach,1,240));
[dep48xFillEach,dep48yFillEach] = avgFill(cat(2,dep48XEach{:}),cat(2,dep48YEach{:}),2);
dep48XMEach = mean(cat(2,dep48XEach{:}),2);
dep48YMEach = mean(cat(2,dep48YEach{:}),2);
[dep48ES,~] = distES(dep48A,reshape(dep48AEach,1,180));
[chowxFillEach,chowyFillEach] = avgFill(cat(2,chowXEach{:}),cat(2,chowYEach{:}),2);
chowXMEach = mean(cat(2,chowXEach{:}),2);
chowYMEach = mean(cat(2,chowYEach{:}),2);
[chowES,~] = distES(chowA,reshape(chowAEach,1,180));
%% Plot generalized and each together
figure
hold on
fill(dep24xFillEach,dep24yFillEach,[0.8 0.8 0.8])
h1 = plot(dep24XMEach,dep24YMEach,'--k','LineWidth',2);
fill(dep24Xfill,dep24Yfill,[0.8 0.8 0.8])
h2 = plot(dep24XM,dep24YM,'-k','LineWidth',2);
legend([h1,h2],{['Individualized: ',num2str(round(mean(reshape(dep24AEach,1,240)),2))],['Generalized: ',num2str(round(mean(dep24A),2))]},'Location','southeast')
text(0.8,0.2,['d = ',num2str(round(dep24ES,2))])
ylim([0 1])
xlim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Dep24: Generalized vs. Individualized')
xlabel('False Positive Rate')
ylabel('True Positive Rate')

figure
hold on
fill(dep48xFillEach,dep48yFillEach,[0.8 0.8 0.8])
h1 = plot(dep48XMEach,dep48YMEach,'--k','LineWidth',2);
fill(dep48Xfill,dep48Yfill,[0.8 0.8 0.8])
h2 = plot(dep48XM,dep48YM,'-k','LineWidth',2);
legend([h1,h2],{['Individualized: ',num2str(round(mean(reshape(dep48AEach,1,180)),2))],['Generalized: ',num2str(round(mean(dep48A),2))]},'Location','southeast')
text(0.8,0.2,['d = ',num2str(round(dep48ES,2))])
ylim([0 1])
xlim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Dep48: Generalized vs. Individualized')
xlabel('False Positive Rate')
ylabel('True Positive Rate')

figure
hold on
fill(chowxFillEach,chowyFillEach,[0.8 0.8 0.8])
h1 = plot(chowXMEach,chowYMEach,'--k','LineWidth',2);
fill(chowXfill,chowYfill,[0.8 0.8 0.8])
h2 = plot(chowXM,chowYM,'-k','LineWidth',2);
legend([h1,h2],{['Individualized: ',num2str(round(mean(reshape(chowAEach,1,180)),2))],['Generalized: ',num2str(round(mean(chowA),2))]},'Location','southeast')
text(0.8,0.2,['d = ',num2str(round(chowES,2))])
ylim([0 1])
xlim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Chow: Generalized vs. Individualized')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
%% Prebinge data
for ii = 1:20
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\preBingeNoise\',num2str(ii),'.mat'])
    preA(ii,:) = concatData{61}.auc;
%     load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\preBingeRand\',num2str(ii),'.mat'])
%     preRandA(ii,:) = concatData{61}.auc;
end
% Plot
figure
hold on
scatterErr(1:61,mean(preA,1),std(preA,[],1),0)
scatterErr(1:61,mean(preRandA,1),std(preRandA,[],1),0)
set(gca,'XTick',1:10:61,'XTickLabel',-2.5:-10:-62.5,'XDir','Reverse','YAxisLocation','right')
xlim([0 62.5])
xlabel('Time before Feeding (sec)')
ylabel('AUC')
title('Predicting Binge Onset')
%% Univariate PreBinge-PreBinge
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\preBingeModelData.mat')
inds = 1:60;
inds = inds(~ismember(inds,[7,17,27,35]));
for ii = 1:20
    disp(ii)
    trainX = allTrainX{ii,1}(:,inds);
    trainY = allTrainY{ii,1};
    for ti = 1:61
        testX = allTestX{ii,ti}(:,inds);
        testY = allTestY{ii,ti};
        for vi = 1:length(inds)
            preMdl{ii,ti,vi} = fitglm(trainX(:,vi),trainY,'distribution','binomial');
            prob = predict(preMdl{ii,ti,vi},testX(:,vi));
            [~,~,~,prePreA(ii,ti,vi)] = perfcurve(testY,prob,1);
        end
    end
end
%% Univariate bingeNot-preBinge data
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\preBingeModelData.mat')
preTestX = allTestX;
preTestY = allTestY;
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\concat50-50_500Train.mat')
%%
inds = 1:60;
inds = inds(~ismember(inds,[7,17,27,35]));
for ii = 1:20
    disp(ii)
    trainX = allTrainX{ii}(:,inds);
    trainY = allTrainY{ii};
    for ti = 1:61
        thisTest = preTestX{ii,ti}(:,inds);
        for vi = 1:length(inds)
            mdl{ii,ti,vi} = fitglm(trainX(:,vi),trainY,'distribution','binomial');
            prob = predict(mdl{ii,ti,vi},thisTest(:,vi));
            [~,~,~,bingeNotPreA(ii,ti,vi)] = perfcurve(preTestY{ii,ti},prob,1);
        end
    end
end
%% Compare univariates
load('C:\USers\Pythia\Documents\GreenLab\data\paper2\analyzed\preVsBinge.mat')
sPrePreA = squeeze(std(prePreA,[],1));
mPrePreA = squeeze(mean(prePreA,1));
[smPrePreA,ord] = sort(mPrePreA(1,:),'descend');
ssPrePreA = sPrePreA(1,ord);

sBingeNotPreA = squeeze(std(bingeNotPreA,[],1));
mBingeNotPreA = squeeze(mean(bingeNotPreA,1));
% [smBingeNotPreA,ord] = sort(mBingeNotPreA,'descend');
ssBingeNotPreA = sBingeNotPreA(ord);

scatterErr(1:56,smPrePreA(1,:),ssPrePreA(1,:),1)
hold on
scatterErr(1:56,mBingeNotPreA(1,ord),sBingeNotPreA(1,ord),0,[0.5 0.5 0.5])
%% Compare beta values from bingeNot and preBinge
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50')
concatBeta = [];
survBeta = [];
for ii = 1:20
    load([num2str(ii),'.mat'])
    concatBeta = [concatBeta;concatData.allBeta{1,1}.betas];
    survBeta = [survBeta;concatData.allBeta{1,1}.survBeta];
end
concatSurv = mean(concatBeta~=0);

% cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\each')
% eachBeta = [];
% survEachBeta = [];
% for ii = 1:240
%     load([num2str(ii),'.mat'])
%     eachBeta = [eachBeta;selfData.allBeta{1,1}.betas];
%     survEachBeta = [survEachBeta;selfData.allBeta{1,1}.survBeta];
% end
% eachSurv = mean(concatBeta~=0);

% Add knonwn NaNs representing values taken out due to noise - see pInds in
% bingeNotData.mat
% concatSurv = [concatSurv(1:6),NaN,concatSurv(7:16),NaN,concatSurv(17:26),NaN,concatSurv(27:34),NaN,concatSurv(35:end)];
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\preBingeNoise')
concatBeta = [];
for ii = 1:20
    load([num2str(ii),'.mat'])
    concatBeta = [concatBeta;concatData{1,1}.allBeta{1,1}.betas];
end
preSurv = mean(concatBeta~=0);
% Stack betas
betas = [concatSurv;preSurv];
%% Apply bingeNot model to preBinge
load('preBingeModelData.mat','allTestX','allTestY','allTestYRand')
inds = 1:60;
inds = inds(~ismember(inds,[7,17,27,35]));
for ii = 1:20
    cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50')
    load([num2str(ii),'.mat'])
    for jj = 1:61
        prob = cvglmnetPredict(concatData.model,zscore(allTestX{ii,jj}(:,inds)),'lambda_1se','response');
        [cpX,cpY,~,cpA(ii,jj)] = perfcurve(allTestY{ii,jj},prob,1);
    end
    cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\500Train_50-50Rand')
    load([num2str(ii),'.mat'])
    for jj = 1:61
        prob = cvglmnetPredict(concatData.model,zscore(allTestX{ii,jj}(:,inds)),'lambda_1se','response');
        [cprX,cprY,~,cprA(ii,jj)] = perfcurve(allTestYRand{ii,jj},prob,1);
    end
end
%%
figure
hold on
scatterErr(1:61,mean(cpA,1),std(cpA,[],1),0)
scatterErr(1:61,mean(cprA,1),std(cprA,[],1),0)
set(gca,'XTick',1:10:61,'XTickLabel',-2.5:-10:-62.5,'XDir','Reverse','YAxisLocation','right')
xlim([0 62.5])
title('BingeNot to PreBinge')
xlabel('Time before Feeding (sec)')
ylabel('AUC')
%%  Apply preBinge model to bingeNot
load('concat50-50_500TrainNew.mat')
for ii = 1:20
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\preBingeNoise\',num2str(ii),'.mat'])
    prob = cvglmnetPredict(concatData{1}.model,zscore(allTestX{ii}(:,inds)),'lambda_1se','response');
    [pbX{ii},pbY{ii},~,pbA(ii)] = perfcurve(allTestY{ii},prob,1);
end
pbRocX = cat(2,pbX{:});
pbRocY = cat(2,pbY{:});
[pbxFill,pbyFill] = avgFill(pbRocX,pbRocY,2);
load('concat50-50_500TrainNewRand.mat')
for ii = 1:20
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\preBingeNoise\',num2str(ii),'.mat'])
    prob = cvglmnetPredict(concatData{1}.model,zscore(allTestX{ii}(:,inds)),'lambda_1se','response');
    [pbXRand{ii},pbYRand{ii},~,pbARand(ii)] = perfcurve(allTestY{ii},prob,1);
end
pbRocXRand = cat(2,pbXRand{:});
pbRocYRand = cat(2,pbYRand{:});
[pbxFillRand,pbyFillRand] = avgFill(pbRocXRand,pbRocYRand,2);
pbES = distES(pbA,pbARand);
% Plot
figure
hold on
fill(pbxFill,pbyFill,[0.8 0.8 0.8])
h(1) = plot(mean(pbRocX,2),mean(pbRocY,2),'k');
fill(pbxFillRand,pbyFillRand,[0.8 0.8 0.8])
h(2) = plot(mean(pbRocXRand,2),mean(pbRocYRand,2),'--k');
h(3) = plot(NaN,NaN,'Color','none');
xlim([0 1])
ylim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
legend([h],{['Real: ',num2str(round(mean(pbA),2))],['Permuted: ',num2str(round(mean(pbARand),2))],['d = ',num2str(round(pbES,2))]},'Location','southeast')
title('PreBinge Model to BingeNot Data')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
%% Power or Coh: PreBinge-Binge-PostBinge
clear chan pair freq
chan = [3];
pair = [];
freq = 3;
feat = 'Low Gamma';
loc = 'SL';
% Get preBinge and notBinge data
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','base','in');
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
   for t = 1:61
       if isempty(pair)
           prePow(ii,t) = mean(psdTrls{t}.relPow(freq,chan,:),'omitnan');
       else
           preCoh(ii,t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
       end
   end
   if isempty(pair)
       notPow(ii) = mean(psdTrls{1,62}.relPow(freq,chan,:),'omitnan');
   else
       notCoh(ii) = mean(coh{1,62}.rel(pair,freq,:));
   end
end

if isempty(pair)
    mPre = mean(prePow,1,'omitnan');
    sPre = std(prePow,[],1,'omitnan');
    mNot = mean(mean(notPow,1,'omitnan'));
else
    mPre = mean(preCoh,1,'omitnan');
    sPre = std(preCoh,[],1,'omitnan');
    mNot = mean(mean(notCoh,1,'omitnan'));
end
% Get Binge data
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge','base','in');
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
    for t = 1:31
        if isempty(pair)
            bingePow(ii,t) = mean(psdTrls{t}.relPow(freq,chan,:),'omitnan');
        else
            bingeCoh(ii,t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
        end
    end
end
if isempty(pair)
    mBinge = mean(bingePow,1,'omitnan');
    sBinge = std(bingePow,[],1,'omitnan');
else
    mBinge = mean(bingeCoh,1,'omitnan');
    sBinge = std(bingeCoh,[],1,'omitnan');
end

% Get PostBinge data
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\postBinge','(e','ex','test','ex');
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
    for t = 1:61
        if isempty(pair)
            postPow(ii,t) = mean(psdTrls{t}.relPow(freq,chan,:),'omitnan');
        else
            postCoh(ii,t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
        end
    end
end
if isempty(pair)
    mPost = mean(postPow,1,'omitnan');
    sPost = std(postPow,[],1,'omitnan');
else
    mPost = mean(postCoh,1,'omitnan');
    sPost = std(postCoh,[],1,'omitnan');
end

% Plot
figure
hold on
shadedErrorBar(1:61,fliplr(mPre.*100),fliplr(sPre.*100))
shadedErrorBar(62:92,fliplr(mBinge.*100),fliplr(sBinge.*100),{'color',[0 0.45 0.74]})
shadedErrorBar(102:110,mPost(1:9).*100,sPost(1:9).*100,{'color',[0 0.45 0.74]})
shadedErrorBar(111:162,mPost(10:61).*100,sPost(10:61).*100)
plot(1:92,ones(1,92).*mNot*100,'k')
plot(1:92,ones(1,92).*mean(mBinge)*100,'--','color',[0 0.45 0.74])
plot(102:162,ones(1,61).*mNot*100,'k')
plot(102:162,ones(1,61).*mean(mBinge)*100,'--','color',[0 0.45 0.74])
xlim([1 162])
set(gca,'XTick',[1:10:51,61.5,71:10:91,101,110.5,121:10:162],'XTickLabel',[-62.5:10:-12.5,0,12.5:10:32.5,-12.5,0,12.5:10:52.5])
title([loc,' ',feat]);
ylabel(['% ',feat])
text(162,mean(mBinge)*100,'Binge','color',[0 0.45 0.74])
text(162,mNot*100,'Other')
xlabel('Time')
box off
%% Core Alpha Power - PreBinge
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','base','in');
for ii = 1:length(files)
   load(files{ii},'psdTrls')
   for t = 1:61
       prePow(ii,t) = mean(psdTrls{t}.relPow(3,4,:));
   end
end
mPre = mean(prePow,1,'omitnan');
sPre = std(prePow,[],1,'omitnan');
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge','base','in');
for ii = 1:length(files)
   load(files{ii},'psdTrls')
   for t = 1:31
       bingePow(ii,t) = mean(psdTrls{t}.relPow(3,4,:));
   end
end
mBinge = mean(bingePow,1,'omitnan');
sBinge = std(bingePow,[],1,'omitnan');
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge','base','in');
for ii = 1:length(files)
   load(files{ii},'psdTrls')
   notPow(ii) = mean(psdTrls{1,2}.relPow(3,4,:));
end
mNot = mean(mean(notPow,1,'omitnan'));
figure
shadedErrorBar(1:31,mBinge.*100,sBinge.*100)
% plot(1:31,mBinge.*100);
hold on
shadedErrorBar(32:92,mPre.*100,sPre.*100)
% plot(32:92,mPre.*100)
plot(1:92,ones(1,92).*mNot*100,'k')
plot(1:92,ones(1,92).*mean(mBinge)*100,'--','color',[0 0.45 0.74])
set(gca,'XTick',[1,11,21,32,44:10:94],'XTickLabel',[32.5,22.5,12.5,0,-12.5:-10:-62.5],'XDir','Reverse','YAxisLocation','right')
title('CR Alpha Power')
ylabel('% Alpha')
xlabel('Time')
text(100,2.1425,'Binge','color',[0 0.45 0.74])
text(101,2.0555,'Other')
box off
%% Core Right Alpha Power - PostBinge
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\postBinge','test','in');
for ii = 1:length(files)
   load(files{ii},'psdTrls')
   for t = 1:61
       prePow(ii,t) = mean(psdTrls{t}.relPow(3,4,:),'omitnan');
   end
end
mPost = mean(prePow,1,'omitnan');
sPost = std(prePow,[],1,'omitnan');
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge','base','in');
for ii = 1:length(files)
   load(files{ii},'psdTrls')
   notPow(ii) = mean(psdTrls{1,2}.relPow(3,4,:),'omitnan');
end
mNot = mean(mean(notPow,1,'omitnan'));
figure
hold on
shadedErrorBar(1:61,mPost.*100,sPost.*100)
plot(1:61,ones(1,61).*mNot*100,'k')
plot(1:61,ones(1,61).*mean(mBinge)*100,'--','color',[0 0.45 0.74])
% set(gca,'XTick',[1,11,21,32,44:10:94],'XTickLabel',[32.5,22.5,12.5,0,-12.5:-10:-62.5],'XDir','Reverse','YAxisLocation','right')
title('CR Alpha Power')
ylabel('% Alpha')
text(101,23.39,'Binge','color',[0 0.45 0.74])
text(101,23.518,'Other')
xlabel('Time')
box off
%% Raw Power
% files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','base','in');
% for ii = 1:length(files)
%    load(files{ii})
%    for t = 1:61
%        pow(ii,t) = mean(psdTrls{t}.bandPow(5,4,:));
%    end
% end
% m = mean(pow,1,'omitnan');
% figure
% plot(1:61,m)
% set(gca,'XTick',1:10:61,'XTickLabel',-2.5:-10:-62.5,'XDir','Reverse','YAxisLocation','right')
% title('CR Low Gamma Power')
% ylabel('dB')
% xlabel('Time')
% %% Total Power
% files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeCombined','base','in');
% for ii = 1:length(files)
%    load(files{ii})
%    for t = 1:61
%        pow(ii,t) = mean(psdTrls{t}.totPow(5,4,:));
%    end
% end
% m = mean(pow,1,'omitnan');
% figure
% plot(1:61,m)
% set(gca,'XTick',1:10:61,'XTickLabel',-2.5:-10:-62.5,'XDir','Reverse','YAxisLocation','right')
% title('CR Total Power')
% ylabel('dB')
% xlabel('Time')
%% Tree analysis
inds = 1:60;
inds = inds(~ismember(inds,[7,17,27,35]));
load('concat50-50_500TrainNew.mat')
tree = [];
impGen = [];
for ii = 1:20
    tree = fitrtree(allTrainX{ii}(:,inds),allTrainY{ii});
    impGen(ii,:) = predictorImportance(tree);
end
load('preBingeModelData.mat','allTrainX','allTrainY')
tree = [];
impPre = [];
for ii = 1:20
   tree = fitrtree(allTrainX{ii,1}(:,inds),allTrainY{ii,1});
   impPre(ii,:) = predictorImportance(tree);
end
figure
scatter(imps(1,:),imps(2,:),100,'.k')
set(gca,'xscale','log','yscale','log')
xlabel('Binge Not')
ylabel('PreBinge')
title('Average Predictor Importance')
%% Build and test univariate logistic from baseline
% Concatenated data
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\concat50-50_500TrainNew.mat')
% inds = 1:60;
% inds = inds(~ismember(inds,pInds));
for ii = 1:20
    trainX = allTrainX{1,ii};
    trainY = allTrainY{1,ii};
    testX = allTestX{1,ii};
    testY = allTestY{1,ii};
    for bi = 1:size(trainX,2)
        mdl{ii,bi} = fitglm(trainX(:,bi),trainY,'distribution','binomial');
        prob = predict(mdl{ii,bi},testX(:,bi));
        [~,~,~,catTestA(ii,bi)] = perfcurve(testY,prob,1);
%         for k = 1:11
%             testX = eachTestX{ii
%             prob = predict(mdl{ii,bi},eachTestX{ii,k}(:,bi));
%             [~,~,~,indTestA(ii,bi,k)] = perfcurve(eachTestY{ii,k},prob,1);
%         end
    end
end
%%
mTestA = mean(catTestA,1);
sTestA = std(catTestA,[],1);
[smTestA,inds] = sort(mTestA,'descend');
ssTestA = sTestA(inds);
scatterErr(1:60,smTestA,ssTestA,1)
hold on
% Plot 50 line
plot([0 60],[0.5 0.5],'--k')
% Plot line of generalized model at baseline (0.7966)
plot([0 60],[mean(ccA) mean(ccA)],'-k')
xlabel('Feature')
ylabel('AUC')
title('Average AUC from Univariate Logistic: Concat')
%% Cat data
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\indvGenCat\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    indCatA(ii,:) = indCatData.auc;
    genA(ii,:) = genData.auc;
    genDep24A(ii,:) = genDep24Data.auc;
    genDep48A(ii,:) = genDep48Data.auc;
    genChowA(ii,:) = genChowData.auc;
    genBaseA(ii,:) = genBaseData.auc;
    baseA(ii,:) = baseData.auc;
    dep24A(ii,:) = dep24Data.auc;
    dep48A(ii,:) = dep48Data.auc;
    chowA(ii,:) = chowData.auc;
end
catInd = [mean(indCatA,1);mean(baseA,1);mean(dep24A,1);mean(dep48A,1);mean(chowA,1)];
catGen = [mean(genA,1);mean(genBaseA,1);mean(genDep24A,1);mean(genDep48A,1);mean(genChowA,1)];

indData = cat(3,indCatA,baseA,dep24A,dep48A,chowA);
genData = cat(3,genA,genBaseA,genDep24A,genDep48A,genChowA);
% for ii = 1:5
%     for jj = 1:8
%         d(ii,jj) = distES(indData(:,jj,ii),genData(:,jj,ii));
%         [~,p(ii,jj)] = ttest(indData(:,jj,ii),genData(:,jj,ii));
%     end
% end
% pNaN = p;
% pNaN(pNaN>0.05) = NaN;
% dNaN = d;
% dNaN(isnan(pNaN)) = NaN;
%%
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\final\indvGenBase\')
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\base500cat8.mat')
baseTestX = eachTestX(:,[1:3,5:7,10,11]);
baseTestY = eachTestY(:,[1:3,5:7,10,11]);
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\pInds.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
for ii = 1:20
    load([num2str(ii),'.mat'])
    indCatA(ii,:) = indCatData.auc;
    genA(ii,:) = genData.auc;
    genDep24A(ii,:) = genDep24Data.auc;
    genDep48A(ii,:) = genDep48Data.auc;
    genChowA(ii,:) = genChowData.auc;
    for jj = 1:8
        [predY] = cvglmnetPredict(baseData.model,zscore(baseTestX{ii,jj}(:,inds)),['lambda_',baseData.hist.cfg.minTerm],'response');
        [~,~,~,genBaseA(ii,jj)] = perfcurve(baseTestY{ii,jj},predY,1);
%         genBaseA(ii,:) = genBaseData.auc;
    end
    baseA(ii,:) = baseData.auc;
    dep24A(ii,:) = dep24Data.auc;
    dep48A(ii,:) = dep48Data.auc;
    chowA(ii,:) = chowData.auc;
end
baseInd = [mean(indCatA,1);mean(baseA,1);mean(dep24A,1);mean(dep48A,1);mean(chowA,1)];
baseGen = [mean(genA,1);mean(genBaseA,1);mean(genDep24A,1);mean(genDep48A,1);mean(genChowA,1)];

indData = cat(3,indCatA,baseA,dep24A,dep48A,chowA);
genData = cat(3,genA,genBaseA,genDep24A,genDep48A,genChowA);
% for ii = 1:5
%     for jj = 1:8
%         d(ii,jj) = distES(indData(:,jj,ii),genData(:,jj,ii));
%         [~,p(ii,jj)] = ttest(indData(:,jj,ii),genData(:,jj,ii));
%     end
% end
% pNaN = p;
% pNaN(pNaN>0.05) = NaN;
% dNaN = d;
% dNaN(isnan(pNaN)) = NaN;
%%
data = [baseInd(1,:)',baseGen(1,:)',catInd(1,:)',catGen(1,:)',baseInd(2,:)',baseGen(2,:)',catInd(2,:)',catGen(2,:)'];
group = {[1,2],[3,4],[1,3],[2,4],[5,6],[7,8],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]};
for ii = 1:size(group,2)
    [~,p(ii)] = ttest(data(:,group{ii}(1)),data(:,group{ii}(2)));
end
pAdj = p*size(p,2);
%%
figure
scatterErr(1:4,mean(data(:,1:4),1),std(data(:,1:4),[],1),0,'k');
scatterErr(5:8,mean(data(:,5:8),1),std(data(:,5:8),[],1),0,[0.5 0.5 0.5]);
sigstar(group(pAdj<0.05),pAdj(pAdj<0.05))
xlim([0.5 8.5])
set(gca,'XTickLabel',{'baseInd','baseGen','catInd','catGen','baseInd','baseGen','catInd','catGen'},'YTick',[0.6:0.1:1])
xlabel('Training Dataset')
ylabel('AUC')
%%
% Plot All states gen vs. ind with averages across animals
figure
h = barwitherr([std(ind,[],2),std(gen,[],2)],[mean(ind,2),mean(gen,2)]);
ylim([0 1.1])
title('Base Data: Gen vs. Ind All States')
set(gca,'xticklabel',{'All','Base','Dep24','Dep48','Chow'})
set(h(1),'facecolor',[0.5 0.5 0.5])
set(h(2),'facecolor','w')
legend({'Ind','Gen'})
ylabel('AUC')
box off
% Plot each state independently across animals (averages are within animal)
% All
figure
h = barwitherr([std(indCatA,[],1);std(genA,[],1)]',[mean(indCatA,1);mean(genA,1)]');
ylim([0 1.1])
title('Base Data: Gen vs. Ind All')
set(h(1),'facecolor',[0.5 0.5 0.5])
set(h(2),'facecolor','w')
legend({'Ind','Gen'})
ylabel('AUC')
xlabel('Animal')
box off
% Base
figure
h = barwitherr([std(baseA,[],1);std(genBaseA,[],1)]',[mean(baseA,1);mean(genBaseA,1)]');
ylim([0 1.1])
title('Base Data: Gen vs. Ind Base')
set(h(1),'facecolor',[0.5 0.5 0.5])
set(h(2),'facecolor','w')
legend({'Ind','Gen'})
ylabel('AUC')
xlabel('Animal')
box off
% Dep24
figure
h = barwitherr([std(dep24A,[],1);std(genDep24A,[],1)]',[mean(dep24A,1);mean(genDep24A,1)]');
ylim([0 1.1])
title('Base Data: Gen vs. Ind Dep24')
set(h(1),'facecolor',[0.5 0.5 0.5])
set(h(2),'facecolor','w')
legend({'Ind','Gen'})
ylabel('AUC')
xlabel('Animal')
box off
% Dep48
figure
h = barwitherr([std(dep48A,[],1);std(genDep48A,[],1)]',[mean(dep48A,1);mean(genDep48A,1)]');
ylim([0 1.1])
title('Base Data: Gen vs. Ind Dep48')
set(h(1),'facecolor',[0.5 0.5 0.5])
set(h(2),'facecolor','w')
legend({'Ind','Gen'})
ylabel('AUC')
xlabel('Animal')
box off
% Chow
figure
h = barwitherr([std(chowA,[],1);std(genChowA,[],1)]',[mean(chowA,1);mean(genChowA,1)]');
ylim([0 1.1])
title('Base Data: Gen vs. Ind Chow')
set(h(1),'facecolor',[0.5 0.5 0.5])
set(h(2),'facecolor','w')
legend({'Ind','Gen'})
ylabel('AUC')
xlabel('Animal')
box off
% Plot p values and Cohen's d
figure
pcolor(padarray(pNaN,[1,1],'post')')
colormap viridis
set(gca,'xtick',1.5:5.5,'xticklabel',{'All','Base','Dep24','Dep48','Chow'},'ytick',1.5:8.5,'yticklabel',1:8)
title('Significant p-values')
figure
pcolor(padarray(dNaN,[1,1],'post')')
colormap viridis
set(gca,'xtick',1.5:5.5,'xticklabel',{'All','Base','Dep24','Dep48','Chow'},'ytick',1.5:8.5,'yticklabel',1:8)
title('Significant Effect Sizes')
%% logCmbs
for ii = 1:3
    disp(num2str(ii))
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\logCmbsFull3\',num2str(ii),'.mat'])
    baseAuc{ii} = A;
    dep24Auc{ii} = dep24A;
    dep48Auc{ii} = dep48A;
    chowAuc{ii} = chowA;
end
%%
baseM = cellfun(@(x) mean(x,1),baseAuc,'UniformOutput',0);
dep24M = cellfun(@(x) mean(x,1),dep24Auc,'UniformOutput',0);
dep48M = cellfun(@(x) mean(x,1),dep48Auc,'UniformOutput',0);
chowM = cellfun(@(x) mean(x,1),chowAuc,'UniformOutput',0);
%% Get top 5 models for each condition and cmbs
for ii  = 1:3
    cmbs{ii} = nchoosek(1:60,ii);
    [sortX,sortInd] = sort(baseM{1,ii},'descend');
    topA(:,ii) = sortX(1:60);
    for jj = 1:60
        topInd{jj,ii} = cmbs{ii}(sortInd(jj),:);
    end
end
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
% inds = 1:60;
% inds = inds(~ismember(inds,[7,17,27,35]));
% nameVect = nameVect(inds);
for ii = 1:60
    for jj = 1:2
        dyad{ii,jj} = nameVect{topInd{ii,2}(jj)};
    end
    for jj = 1:3
        triad{ii,jj} = nameVect{topInd{ii,3}(jj)};
    end
end
%%
[~,ord] = sort(baseM{3},'descend');
figure
subplot(2,2,1)
hold on
scatter(1:length(baseM{3}),sort(baseM{3},'descend'),'.')
scatter(1:length(dep24M{3}),sort(dep24M{3},'descend'),'.')
scatter(1:length(dep48M{3}),sort(dep48M{3},'descend'),'.')
scatter(1:length(chowM{3}),sort(chowM{3},'descend'),'.')

subplot(2,2,2)
hold on
scatter(1:34220,baseM{3}(ord),'.')
scatter(1:34220,dep24M{3}(ord),1,[0.85 0.325 0.098])
scatter(1:34220,sort(dep24M{3},'descend'),'.k')

subplot(2,2,3)
hold on
scatter(1:34220,baseM{3}(ord),'.')
scatter(1:34220,dep48M{3}(ord),1,[0.929 0.694 0.125])
scatter(1:34220,sort(dep48M{3},'descend'),'.k')

subplot(2,2,4)
hold on
scatter(1:34220,baseM{3}(ord),'.')
scatter(1:34220,chowM{3}(ord),1,[0.494 0.184 0.556])
scatter(1:34220,sort(chowM{3},'descend'),'.k')
%%
baseAM = cell2mat(cellfun(@(x) mean(mean(x,1)),baseAuc,'UniformOutput',0));
baseAS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),baseAuc,'UniformOutput',0));
dep24AM = cell2mat(cellfun(@(x) mean(mean(x,1)),dep24Auc,'UniformOutput',0));
dep24AS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),dep24Auc,'UniformOutput',0));
dep48AM = cell2mat(cellfun(@(x) mean(mean(x,1)),dep48Auc,'UniformOutput',0));
dep48AS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),dep48Auc,'UniformOutput',0));
chowAM = cell2mat(cellfun(@(x) mean(mean(x,1)),chowAuc,'UniformOutput',0));
chowAS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),chowAuc,'UniformOutput',0));
%% Plot
figure
subplot(2,2,1)
scatterErr(1:58,baseAM,baseAS,0)
ylim([0.5 0.8])
title('Baseline-Baseline Model')
subplot(2,2,2)
scatterErr(1:58,dep24AM,dep24AS,0)
ylim([0.5 0.8])
title('Baseline-Dep24 Model')
subplot(2,2,3)
scatterErr(1:58,dep48AM,dep48AS,0)
ylim([0.5 0.8])
title('Baseline-Dep48 Model')
subplot(2,2,4)
scatterErr(1:58,chowAM,chowAS,0)
ylim([0.5 0.8])
title('Baseline-Chow Model')
%% Test best model across times
best = inds(2);
for ii = 1:size(allTestX,1)
   for jj = 1:size(allTestX,2)
       prob = predict(mdl{ii,best},allTestX{ii,jj}(:,best));
       [~,~,~,a(ii,jj)] = perfcurve(allTestY{ii,jj},prob,1);
   end
end
%%
figure
scatterErr(0:60,mean(a,1),std(a,[],1),1)
set(gca,'XTickLabel',2.5:10:62.5)
xlabel('Time before Eating (sec)')
ylabel('AUC')
%% 
load('I6Base_2015-11-24.mat')
disp('Applying 60 Hz filter with filter60.m...')
[LFPTs.data] = filter60(LFPTs,adfreq,'off');
disp('Downsampling signal with dwnSample.m...')
[LFPTs,adfreq] = dwnSample(LFPTs,5,adfreq);
%%
% Calculate spectrogram
[s,f,t] = spectrogram(LFPTs.data(1,:),2048,1024,1:100,400);
% Convert s to dB
pow = 10*log10(abs(s));
% For each window (column of pow) get and normalize high gamma power
hgPowNorm = trapz(pow(70:90,:),1)./trapz(pow(1:100,1));
% Plot normalized high gamma power
figure
hold on
plot(t,hgPowNorm)
% Add binge events
for ii = 1:size(eventTs.t{1,7},1)
    plot([eventTs.t{1,7}(ii) eventTs.t{1,8}(ii)],[1 1],'-k','LineWidth',2)
end
%% Extract all shell left high gamma (ind = 6) from all test sets
preOne = cell(1,61); preZero = cell(1,61);
for ii = 1:size(allTestY,1)
    for jj = 1:size(allTestY,2)
        % Get inds of ones
        oneInd{ii,jj} = logicFind(1,allTestY{ii,jj},'==');
        % Get inds of zeros
        zeroInd{ii,jj} = logicFind(0,allTestY{ii,jj},'==');
        % Use above to extract and concatenate shell left high gamma
        preOne{jj} = [preOne{jj};allTestX{ii,jj}(oneInd{ii,jj},6)];
        preZero{jj} = [preZero{jj};allTestX{ii,jj}(zeroInd{ii,jj},6)];
    end
end
%% Average each time step
preAvg(1,:) = cellfun(@mean,preOne,'UniformOutput',0);
preAvg(2,:) = cellfun(@mean,preZero,'UniformOutput',0);
%%
% Extract binge start time using second column of samp
ind = cell(1);
for k = 1:size(samp{1,1},1)
    starts = [];
    for ii = 1:61
        starts = [starts;samp{1,1}{k,ii}(:,2)+(400*(ii-1))];
    end
    starts = unique(starts);
    for jj = 1:length(starts)
        for ii = 1:61
            thisInd = logicFind(starts(jj)-400*(ii-1),samp{1,1}{k,ii}(:,2),'==');
            if isempty(thisInd)
                ind{k}(jj,ii) = NaN;
            else
                ind{k}(jj,ii) = thisInd;
            end
        end
    end
end
%% Use indices to extract shell left high gamma from pre-trials
for ii = 1:23
    for jj = 1:61
        for k = 1:size(ind{ii},1)
            if isnan(ind{ii}(k,jj))
                hgPow{ii}(k,jj) = NaN;
            else
                hgPow{ii}(k,jj) = data{1,1}{ii,jj}(ind{ii}(k,jj),6);
            end
        end
    end
    hgPowRand{ii} = data{1,1}{ii,end}(:,6);
end
%%
test = diff(samp{1,1}{1,62}(:,2));
contig = [];
contig(1) = 1;
r = 2;
c = 1;
ii = 2;
while ii < size(samp{1,1}{1,62},1)
    if samp{1,1}{1,62}(ii,1) == samp{1,1}{1,62}(ii-1,1)+2000
        c = c+1;
        contig(r-1,c) = ii;
        ii = ii+1;
    else
        c = 1;
        r = r+1;
        contig(r-1,c) = ii;
        ii = ii+1;
    end
end
%%
allPrePow = cat(1,hgPow{:});
allRandPow = 
%%
% Get average ROC for concat-concat models
ccMX = mean(cat(2,ccX{:}),2);
ccMY = mean(cat(2,ccY{:}),2);
% Get average 'fill' for concat-concat models
[ccXfill,ccYfill] = avgFill(cat(2,ccX{:}),cat(2,ccY{:}),2,1);
%% Combine other conditions and plot
for ii = 1:57
   m(ii) = mean([reshape(dep24Auc{1,ii},1,numel(dep24Auc{1,ii})),reshape(dep48Auc{1,ii},1,numel(dep48Auc{1,ii})),reshape(chowAuc{1,ii},1,numel(chowAuc{1,ii}))]); 
   s(ii) = std([reshape(dep24Auc{1,ii},1,numel(dep24Auc{1,ii})),reshape(dep48Auc{1,ii},1,numel(dep48Auc{1,ii})),reshape(chowAuc{1,ii},1,numel(chowAuc{1,ii}))]); 
end
scatterErr(1:57,m,s,1)
%% Load and analyze given animals/binge epochs
fNames = {'I2BaseDec15.mat','I1BaseNov9.mat','I12BaseNov12.mat','H15BaseSep26.mat'};
trial = [1,35,5,8];
for ii = 1:4
    load(fNames{ii})
    start{ii} = nearest_idx3(eventTs.t{1,7}(trial(ii)),LFPTs.tvec);
    stop{ii} = nearest_idx3(eventTs.t{1,8}(trial(ii)),LFPTs.tvec);
    x1 = LFPTs.data(4,start{ii}-20000:start{ii}+20000);
    y1 = LFPTs.data(3,start{ii}-20000:start{ii}+20000);
    x1 = filter60(x1,adfreq,'off');
    y1 = filter60(y1,adfreq,'off');
    [pow{ii},f{ii},coi{ii}] = cwt(x1,'morse',adfreq);
    [wcoh{ii},~,cohf{ii}] = wcoherence(x1,y1,adfreq);
    %     [s{ii},f{ii},~] = spectrogram(LFPTs.data(4,start{ii}-20000:start{ii}+20000),512,256,1:200,adfreq);
    x2 = LFPTs.data(4,start{ii}-4000:start{ii}+4000);
    y2 = LFPTs.data(3,start{ii}-4000:start{ii}+4000);
    x2 = filter60(x2,adfreq,'off');
    y2 = filter60(y2,adfreq,'off');
    [pow2{ii},f2{ii},coi2{ii}] = cwt(x2,'morse',adfreq);
    [wcoh2{ii},~,cohf2{ii}] = wcoherence(x2,y2,adfreq);
    %     [s2{ii},f2{ii},~] = spectrogram(LFPTs.data(4,start{ii}-4000:start{ii}+4000),256,128,1:200,adfreq);
    t{ii} = (start{ii}-20000:start{ii}+20000)./adfreq;
    t2{ii} = (start{ii}-4000:start{ii}+4000)./adfreq;
end
%% Plot wavelet PSD
c = [0.06,0.08,0.07,0.24];
for ii = 1:4
    figure
    ftest = log2(f{ii});
%     h = pcolor(t{ii},ftest,abs(pow{ii}));
%     h.EdgeColor = 'none';
    imagesc(t{ii},ftest,abs(pow{ii}));
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    caxis([0 c(ii)])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    
    title([fNames{ii}])
    
    figure
    ftest = log2(f2{ii});
%     h = pcolor(t2{ii},ftest,abs(pow2{ii}));
%     h.EdgeColor = 'none';
    imagesc(t2{ii},ftest,abs(pow2{ii}));
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    caxis([0 c(ii)])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    title([fNames{ii}])
end
%% Plot wavelet COH
for ii = 1:4
    figure
    ftest = log2(cohf{ii});
%     h = pcolor(t{ii},ftest,wcoh{ii});
%     h.EdgeColor = 'none';
    imagesc(t{ii},ftest,wcoh{ii})
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    title([fNames{ii}])
    
    figure
    ftest = log2(cohf2{ii});
%     h = pcolor(t2{ii},ftest,wcoh2{ii});
%     h.EdgeColor = 'none';
    imagesc(t2{ii},ftest,wcoh2{ii})
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    title([fNames{ii}])
end
%% Subset Indices
% Power indices (channel,band)
powInd = reshape(1:24,6,4)';
% Coherence indices (cmb,band)
cohInd = reshape(25:60,6,6)';
% Delta indices
dInd = [powInd(:,1);cohInd(:,1)]';
% Theta indices
tInd = [powInd(:,2);cohInd(:,2)]';
% Alpha indices
aInd = [powInd(:,3);cohInd(:,3)]';
% Beta indices
bInd = [powInd(:,4);cohInd(:,4)]';
% Low gamma indices
lgInd = [powInd(:,5);cohInd(:,5)]';
% High gamma indices
hgInd = [powInd(:,6);cohInd(:,6)]';
% Left indices
lInd = [powInd([1,3],:);cohInd(2,:)];
% Right indices
rInd = [powInd([2,4],:);cohInd(5,:)];
% Shell indices
sInd = [powInd([1,2],:);cohInd(1,:)];
% Core indices
cInd = [powInd([3,4],:);cohInd(6,:)];
%
inds{1} = powInd;
inds{2} = cohInd;
inds{3} = dInd;
inds{4} = tInd;
inds{5} = aInd;
inds{6} = bInd;
inds{7} = lgInd;
inds{8} = hgInd;
inds{9} = lInd;
inds{10} = rInd;
inds{11} = sInd;
inds{12} = cInd;
%% Load subset files
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\concat\subset\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    for jj = 1:12
       if isequal(concatData.rocX{jj},[0;1])
          concatData.rocX{jj} = (0:1/1500:1)';
          concatData.rocY{jj} = (0:1/1500:1)';
       end
    end
    ccSX(:,:,ii) = cat(2,concatData.rocX{:});
    ccSY(:,:,ii) = cat(2,concatData.rocY{:});
    ccSA(ii,:) = cell2mat(concatData.auc);
end
% Load subsetRand files
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\concat\subsetRand\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    for jj = 1:12
       if isequal(concatData.rocX{jj},[0;1])
          concatData.rocX{jj} = (0:1/1500:1)';
          concatData.rocY{jj} = (0:1/1500:1)';
       end
    end
    ccSRandX(:,:,ii) = cat(2,concatData.rocX{:});
    ccSRandY(:,:,ii) = cat(2,concatData.rocY{:});
    ccSRandA(ii,:) = cell2mat(concatData.auc);
end
% Get average fills
for ii = 1:12
    [ccSXfill(ii,:),ccSYfill(ii,:)] = avgFill(ccSX(:,ii,:),ccSY(:,ii,:),3);
    [ccSXrandFill(ii,:),ccSYrandFill(ii,:)] = avgFill(ccSRandX(:,ii,:),ccSRandY(:,ii,:),3);
end
% Get average subset AUC
mccSA = mean(ccSA,1);
sccSA = std(ccSA,[],1);
%% Plot each subset
figure
hold on
for ii = 1:12
   plot(mean(ccSX(:,ii,:),3),mean(ccSY(:,ii,:),3)) 
end
% Add concat-concat line
plot(ccMX,ccMY,'-k','LineWidth',2)
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Subset to Concat ROCs')
%%
figure
for ii = 1:12
    subplot(4,3,ii)
    hold on
    fill(ccXfill(ii,:),ccYfill(ii,:),'b')
    fill(ccXrandFill(ii,:),ccYrandFill(ii,:),'b')
    xlim([0 1])
    ylim([0 1])
end
%% Get effect size for each subset
for ii = 1:12
   [subsetES(ii),~] = distES(ccA(:,ii),ccRandA(:,ii)); 
end
