[dataBlast,sampBlast,filesBlast] = collateData('E:\blast\processed\blast\',...
    {'.mat'},{'pow','coh'},'trl','rel');
[dataCon,sampCon,filesCon] = collateData('E:\blast\processed\control\',...
    {'.mat'},{'pow','coh'},'trl','rel');
%%
this = cellfun(@(x) strsplit(x,'_'),filesBlast{1},'UniformOutput',0);
for ii = 1:numel(this)
    blastName(ii) = this{ii}(1);
    blast{ii} = dataBlast{1}{ii};
end
uBlast = unique(blastName);

this = cellfun(@(x) strsplit(x,'_'),filesCon{1},'UniformOutput',0);
for ii = 1:numel(this)
    conName(ii) = this{ii}(1);
    con{ii} = dataCon{1}{ii};
end
uCon = unique(conName);
%% low and high gamma plots
lg = [5,11,17,23,29,35,41,47];
hg = [6,12,18,24,30,36,42,48];
thisCon = cat(1,con{:});
thisBlast = cat(1,blast{:});
thisConM = cellfun(@(x) mean(x,1),con,'UniformOutput',false);
thisConM = cat(1,thisConM{:});
thisBlastM = cellfun(@(x) mean(x,1),blast,'UniformOutput',false);
thisBlastM = cat(1,thisBlastM{:});
for ii = 1:numel(lg)
    [~,lgP(ii)] = ttest2(thisConM(:,lg(ii)),thisBlastM(:,lg(ii)));
    doubleHist(thisConM(:,lg(ii)),thisBlastM(:,lg(ii)))    
end
for ii = 1:numel(hg)
    [~,hgP(ii)] = ttest2(thisConM(:,hg(ii)),thisBlastM(:,hg(ii)));
    doubleHist(thisConM(:,hg(ii)),thisBlastM(:,hg(ii)))    
end
%% LOO
minSamp = 2000;
cmbs = [];
c = 1;
for ii = 1:numel(uCon)
    for jj = 1:numel(uBlast)
        cmbs(c,:) = [ii,jj];
        c = c+1;
    end
end
for ii = 1:size(cmbs,1)
    disp(ii)
    thisBlast = [];
    for k = 1:numel(uBlast)
        blastInd = logicFind(uBlast{k},blastName,'==');
        theseBlast = cat(1,blast{blastInd});
        thisBlast{k} = theseBlast(randperm(size(theseBlast,1),...
            minSamp),:);
    end
    thisCon = [];
    for k = 1:numel(uCon)
        conInd = logicFind(uCon{k},conName,'==');
        theseCon = cat(1,con{conInd});
        thisCon{k} = theseCon(randperm(size(theseCon,1),...
            minSamp),:);
    end

    otherBlastInd = logicFind(1,~ismember(1:numel(uBlast),...
        cmbs(ii,1)),'==');
    leftBlast = cat(1,thisBlast{cmbs(ii,1)});

    otherConInd = logicFind(1,~ismember(1:numel(uCon),...
        cmbs(ii,2)),'==');
    leftCon = cat(1,thisCon{cmbs(ii,2)});

    mdl = fitglm(cat(1,thisBlast{otherBlastInd},...
        thisCon{otherConInd}),...
        cat(1,ones(size(cat(1,thisBlast{otherBlastInd}),1),1),...
        zeros(size(cat(1,thisCon{otherConInd}),1),1)),...
        'distribution','binomial');
    prob = predict(mdl,cat(1,leftBlast,leftCon));
    [~,~,~,aLOO(ii)] = perfcurve(cat(1,ones(size(leftBlast,1),1),...
        zeros(size(leftCon,1),1)),prob,1);
    % Permuted
    blastPermInd = randperm(11,6);
    blastPerm = cat(1,thisBlast{blastPermInd});
    otherBlastInd = logicFind(1,~ismember(1:11,blastPermInd),'==');
    otherBlastPerm = cat(1,thisBlast{otherBlastInd});

    conPermInd = randperm(11,6);
    conPerm = cat(1,thisCon{conPermInd});
    otherConInd = logicFind(1,~ismember(1:11,conPermInd),'==');
    otherConPerm = cat(1,thisCon{otherConInd});

    mdl = fitglm(cat(1,blastPerm,conPerm,otherBlastPerm,otherConPerm),...
        cat(1,ones(size(cat(1,blastPerm,conPerm),1),1),...
        zeros(size(cat(1,otherBlastPerm,otherConPerm),1),1)),...
        'Distribution','binomial');
    prob = predict(mdl,cat(1,leftBlast,leftCon));
    [~,~,~,aLOOP(ii)] = perfcurve(cat(1,ones(size(leftBlast,1),1),...
        zeros(size(leftCon,1),1)),prob,1);
end
%% 80:20
minSamp = 2000;
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
for ii = 1:200
    disp(ii)
    thisBlast = [];
    for k = 1:numel(uBlast)
        blastInd = logicFind(uBlast{k},blastName,'==');
        theseBlast = cat(1,blast{blastInd});
        thisBlast{k} = theseBlast(randperm(size(theseBlast,1),...
            minSamp),:);
    end
    thisCon = [];
    for k = 1:numel(uCon)
        conInd = logicFind(uCon{k},conName,'==');
        theseCon = cat(1,con{conInd});
        thisCon{k} = theseCon(randperm(size(theseCon,1),...
            minSamp),:);
    end
    [trainX,trainY,testX,testY] = trainTest(cat(1,thisBlast{:},...
        thisCon{:}),cat(1,ones(12*minSamp,1),zeros(12*minSamp,1)),0.2);
    mdl = fitglm(trainX,trainY,'Distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,a80(ii)] = perfcurve(testY,prob,1);
    for f = 1:216
        mdl = fitglm(trainX(:,f),trainY,'Distribution','binomial');
        prob = predict(mdl,testX(:,f));
        [~,~,~,a80s(ii,f)] = perfcurve(testY,prob,1);
        coeffs(ii,f) = table2array(mdl.Coefficients(2,1));
    end
    % low gamma permute
    lgInds = logicFind(1,contains(feat(1:216),'lg'),'==');
    newTrainX = trainX;
    for k = lgInds
        newTrainX(:,k) = trainX(randperm(size(trainX,1),size(trainX,1)),k);
    end
    mdl = fitglm(newTrainX,trainY,'Distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,a80LG(ii)] = perfcurve(testY,prob,1);
    % hg permute
    hgInds = logicFind(1,contains(feat(1:216),'hg'),'==');
    newTrainX = trainX;
    for k = hgInds
        newTrainX(:,k) = trainX(randperm(size(trainX,1),size(trainX,1)),k);
    end
    mdl = fitglm(newTrainX,trainY,'Distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,a80HG(ii)] = perfcurve(testY,prob,1);
    % lg and hg permute
    gInds = cat(2,lgInds,hgInds);
    newTrainX = trainX;
    for k = gInds
        newTrainX(:,k) = trainX(randperm(size(trainX,1),size(trainX,1)),k);
    end
    mdl = fitglm(newTrainX,trainY,'Distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,a80G(ii)] = perfcurve(testY,prob,1);
    % Permuted
    blast1 = randperm(12,6);
    blast2 = logicFind(1,~ismember(1:12,blast1),'==');
    con1 = randperm(12,6);
    con2 = logicFind(1,~ismember(1:12,con1),'==');
    
    group1 = cat(1,thisBlast{blast1},thisCon{con1});
    group2 = cat(1,thisBlast{blast2},thisCon{con2});

    [trainX,trainY,testX,testY] = trainTest(cat(1,group1,group2),...
        cat(1,ones(size(group1,1),1),zeros(size(group2,1),1)),0.2);
    mdl = fitglm(trainX,trainY,'Distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,a80p(ii)] = perfcurve(testY,prob,1);
end
% save('E:/blast/blastConRelModels.mat','a80','a80p','a80s','coeffs')
%%
% figure
% subplot(2,1,1)
% [f,xi,bw] = ksdensity(aLOO);
% plot(xi,f.*bw)
% % fill(xi,f/100,'w')
% xlim([0 1])
% title('blast vs. con: LOO')
% 
% subplot(2,1,2)
% [f,xi,bw] = ksdensity(aLOOP);
% hold on
% plot(xi,f.*bw)
% plot([mean(aLOO) mean(aLOO)],[0 .15])
% % fill(xi,f/100,'w')
% xlim([0 1])
% title('blast vs. con: LOO permuted')

figure
subplot(2,1,1)
[f,xi,bw] = ksdensity(a80);
plot(xi,f.*bw)
% fill(xi,f/100,'w')
hold on
% [f,xi,bw] = ksdensity(a80LG);
% plot(xi,f.*bw)
% [f,xi,bw] = ksdensity(a80HG);
% plot(xi,f.*bw)
% [f,xi,bw] = ksdensity(a80G);
% plot(xi,f.*bw)
plot([mean(a80LG) mean(a80LG)],[0 .15])
plot([mean(a80HG) mean(a80HG)],[0 .15])
plot([mean(a80G) mean(a80G)],[0 .15])
xlim([0.5 1])
title('blast vs. con: 80:20')
legend({'full','LG permute','HG permute','G permute'})

subplot(2,1,2)
[f,xi,bw] = ksdensity(a80p);
hold on
plot(xi,f.*bw,'k')
plot([mean(a80) mean(a80)],[0 .15],'color',[0 0.4470 0.7410])
% fill(xi,f/100,'w')
title('blast vs. con: 80:20 permuted')
xlim([0.5 1])
%% single feature
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
[sA80s,ind] = sort(mean(a80s,1),'descend');
featSort = feat(ind)';
figure
plot(1:216,sA80s,'o')
%% correlate top feature with behavioral data
% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 7);
% Specify sheet
opts.Sheet = "bTBI + ephys";
% Specify column names and types
opts.VariableNames = ["Var1", "VarName2", "contextualFearConditioning", "OpenField", "Var5", "SucrosePreference", "MorrisWaterMaze"];
opts.SelectedVariableNames = ["VarName2", "contextualFearConditioning", "OpenField", "SucrosePreference", "MorrisWaterMaze"];
opts.VariableTypes = ["char", "string", "double", "double", "char", "double", "double"];
% Specify variable properties
opts = setvaropts(opts, ["Var1", "VarName2", "Var5"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "VarName2", "Var5"], "EmptyFieldRule", "auto");
% Import the data
behav = table;
ranges = ["A3:G8", "A10:G16"];
for idx = 1:length(ranges)
    opts.DataRange = ranges(idx);
    tb = readtable("H:\Shared drives\dwielDoucetteLab\data\blast\Blast_correlation_ephys+MRI_[A1].xlsx", opts, "UseExcel", false);
    behav = [behav; tb]; %#ok<AGROW>
end
load('blastConRelModels2.mat')
load('blastCon.mat')
%%
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
[sA80s,ind] = sort(mean(a80s,1),'descend');
featSort = feat(ind)';
%%
thisData = zeros(size(behav,1),216);
for ii = 1:size(behav,1)
    if ii <= 6
        theseInds = contains(blastName,behav.VarName2(ii));
        thisData(ii,:) = mean(cat(1,dataBlast{1}{theseInds}),1);
    else
        theseInds = contains(conName,behav.VarName2(ii));
        thisData(ii,:) = mean(cat(1,dataCon{1}{theseInds}),1);
    end
end
%%
for ii = 1:5
    figure
    subplot(2,2,1)
    plot(thisData(:,ind(ii)),behav.contextualFearConditioning,'.')
    lsline
    mdl = fitlm(thisData(:,ind(ii)),behav.contextualFearConditioning);
    logMdl = fitglm(thisData(:,ind(ii)),...
        behav.contextualFearConditioning>=...
        median(behav.contextualFearConditioning),...
        'distribution','binomial');
    title(['fear R^2=',num2str(mdl.Rsquared.Ordinary),' p=',...
        num2str(table2array(logMdl.Coefficients(2,4)))])

    subplot(2,2,2)
    plot(thisData(:,ind(ii)),behav.OpenField,'.')
    lsline
    mdl = fitlm(thisData(:,ind(ii)),behav.OpenField);
    logMdl = fitglm(thisData(:,ind(ii)),...
        behav.OpenField>=...
        median(behav.OpenField),...
        'distribution','binomial');
    title(['open field R^2=',num2str(mdl.Rsquared.Ordinary),' p=',...
        num2str(table2array(logMdl.Coefficients(2,4)))]);

    subplot(2,2,3)
    plot(thisData(:,ind(ii)),behav.SucrosePreference,'.')
    lsline
    mdl = fitlm(thisData(:,ind(ii)),behav.SucrosePreference);
    logMdl = fitglm(thisData(:,ind(ii)),...
        behav.SucrosePreference>=...
        median(behav.SucrosePreference),...
        'distribution','binomial');
    title(['sucrose R^2=',num2str(mdl.Rsquared.Ordinary),' p=',...
        num2str(table2array(logMdl.Coefficients(2,4)))])

    subplot(2,2,4)
    plot(thisData(:,ind(ii)),behav.MorrisWaterMaze,'.')
    lsline
    mdl = fitlm(thisData(:,ind(ii)),behav.MorrisWaterMaze);
    logMdl = fitglm(thisData(:,ind(ii)),...
        behav.MorrisWaterMaze>=...
        median(behav.MorrisWaterMaze),...
        'distribution','binomial');
    title(['MWM R^2=',num2str(mdl.Rsquared.Ordinary),' p=',...
        num2str(table2array(logMdl.Coefficients(2,4)))])

    sgtitle(feat(ind(ii)))
end
%% correlate top features with max freezing day 1 and total freezing day 2
freeze1 = [70.0002,60.6946,31.9442,32.9168,96.3889,94.444,83.056,100,100,59.7222,96.9444,72.4997];
freeze2 = [0.208,39.0386,1.2601,4.1436,17.8998,7.9676,0.7041,14.7616,11.6484,61.0427,29.5445,12.2613];
thisFeat = ind(1);
for ii = 1:6
    % 1d2, 3d1, etc.
    this{1,ii} = mean(cat(1,dataBlast{1}{4*ii-3:4*ii,:}));
    % VA10, VA16, etc.
    this{2,ii} = mean(cat(1,dataBlast{1}{3*ii-2:3*ii,:}));
end
figure
subplot(2,2,1)
hold on
for ii = 1:6
    that(1,ii) = this{2,ii}(thisFeat);
    that(2,ii) = freeze1(ii);
end
plot(that(1,:),that(2,:),'.k')
lsline
mdl = fitlm(that(1,:),that(2,:));
title(['VA vs. freeze day 1; R^2 = ',num2str(mdl.Rsquared.Ordinary)])

subplot(2,2,2)
hold on
for ii = 1:6
    that(1,ii) = this{2,ii}(thisFeat);
    that(2,ii) = freeze2(ii);
end
plot(that(1,:),that(2,:),'.k')
lsline
mdl = fitlm(that(1,:),that(2,:));
title(['VA vs. freeze day 2; R^2 = ',num2str(mdl.Rsquared.Ordinary)])

subplot(2,2,3)
hold on
for ii = 1:6
    that(1,ii) = this{1,ii}(thisFeat);
    that(2,ii) = freeze1(ii+6);
end
plot(that(1,:),that(2,:),'.k')
lsline
mdl = fitlm(that(1,:),that(2,:));
title(['d vs. freeze day 1; R^2 = ',num2str(mdl.Rsquared.Ordinary)])

subplot(2,2,4)
hold on
for ii = 1:6
    that(1,ii) = this{1,ii}(thisFeat);
    that(2,ii) = freeze2(ii+6);
end
plot(that(1,:),that(2,:),'.k')
lsline
mdl = fitlm(that(1,:),that(2,:));
title(['d vs. freeze day 2; R^2 = ',num2str(mdl.Rsquared.Ordinary)])