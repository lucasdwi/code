[dataBlast,sampBlast,filesBlast] = collateData('E:\blast\processed\blast\',...
    {'.mat'},{'pow','coh'},'trl','rel');
[dataCon,sampCon,filesCon] = collateData('E:\blast\processed\control\',...
    {'.mat'},{'pow','coh'},'trl','rel');
%%
this = cellfun(@(x) strsplit(x,'_'),filesBlast{1},'UniformOutput',0);
for ii = 1:numel(this)
    blastName(ii) = this{ii}(1);
%     blastZ{ii} = zscore(dataBlast{1}{ii});
    blastZ{ii} = dataBlast{1}{ii};
end
uBlast = unique(blastName);

this = cellfun(@(x) strsplit(x,'_'),filesCon{1},'UniformOutput',0);
for ii = 1:numel(this)
    conName(ii) = this{ii}(1);
%     conZ{ii} = zscore(dataCon{1}{ii});
    conZ{ii} = dataCon{1}{ii};
end
uCon = unique(conName);
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
        theseBlast = cat(1,blastZ{blastInd});
        thisBlast{k} = theseBlast(randperm(size(theseBlast,1),...
            minSamp),:);
    end
    thisCon = [];
    for k = 1:numel(uCon)
        conInd = logicFind(uCon{k},conName,'==');
        theseCon = cat(1,conZ{conInd});
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
for ii = 101:200
    disp(ii)
    thisBlast = [];
    for k = 1:numel(uBlast)
        blastInd = logicFind(uBlast{k},blastName,'==');
        theseBlast = cat(1,blastZ{blastInd});
        thisBlast{k} = theseBlast(randperm(size(theseBlast,1),...
            minSamp),:);
    end
    thisCon = [];
    for k = 1:numel(uCon)
        conInd = logicFind(uCon{k},conName,'==');
        theseCon = cat(1,conZ{conInd});
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
save('E:/blast/blastConModels.mat','aLOO','aLOOP','a80','a80p','a80s','coeffs')
%%
figure
subplot(2,1,1)
[f,xi,bw] = ksdensity(aLOO);
plot(xi,f.*bw)
% fill(xi,f/100,'w')
xlim([0 1])
title('blast vs. con: LOO')

subplot(2,1,2)
[f,xi,bw] = ksdensity(aLOOP);
hold on
plot(xi,f.*bw)
plot([mean(aLOO) mean(aLOO)],[0 .15])
% fill(xi,f/100,'w')
xlim([0 1])
title('blast vs. con: LOO permuted')

figure
subplot(2,1,1)
[f,xi,bw] = ksdensity(a80);
plot(xi,f.*bw)
% fill(xi,f/100,'w')
xlim([0 1])
title('blast vs. con: 80:20')

subplot(2,1,2)
[f,xi,bw] = ksdensity(a80p);
hold on
plot(xi,f.*bw)
plot([mean(a80) mean(a80)],[0 .15])
% fill(xi,f/100,'w')
title('blast vs. con: 80:20 permuted')
xlim([0 1])