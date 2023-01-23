% state models
load('F:/irdmRound2/ddtFiles.mat')
files = fileSearch('F:/irdmRound2/processedContinuous3/','.mat');
inter = {'sIL','sNAcC','guan0.3','mph3','mph0.3','mph1','mph0.1',...
    'sNAcC-LSD','sIL-LSD'};
interInds = cellfun(@(x) contains(processedData(:,6),x),...
    inter,'UniformOutput',false);
interInds{1}(interInds{9}) = 0;
interInds{2}(interInds{8}) = 0;
% any(sum(cat(2,interInds{:}),2)>1);
for ii = 1:numel(inter)
    interInds{ii} = interInds{ii}.*ii;
end
interGroup = sum(cat(2,interInds{:}),2);
clear

for ii = 11:size(processedData,1)
    disp([num2str(ii),' of ',num2str(size(processedData,1))])
    load(files{contains(files,processedData{ii,1}) & ...
        contains(files,processedData{ii,2})},'LFPTs','psdTrls','coh');
    pow = squeeze(mean(psdTrls{1}.relPow,4,'omitnan'));
    coh = squeeze(mean(coh{1}.normBandCoh,4,'omitnan')); 
    % Check for NaNs;
    if size(LFPTs.data,1) < 8
        [pow,coh] = addNaNs(LFPTs,pow,coh,8);
    end
    rePow = reshape(pow,48,1)';
    reCoh = reshape(permute(coh,[2,1,3]),168,1)';
    
    data{ii,1} = [rePow,reCoh];
end
data(:,2:4) = processedData(:,[1:2,8]);
data(:,5) = num2cell(interGroup);
% save('F:/irdmRound2/stateData.mat','data')
%% Add in sex info; 1 = male and 0 = female
sexData = readtable('F:/irdmRound2/sex.xlsx');
for ii = 1:size(data,1)
    this = str2num(data{ii,2}(end-1:end));
    sexInd = logicFind(this,table2array(sexData(:,1)),'==');
    if contains(table2cell(sexData(sexInd,2)),'M')
        data{ii,6} = 1;
    else
        data{ii,6} = 0;
    end
end
% save('F:/irdmRound2/stateData.mat','data','sexData')
%% Build base model, LOO and 80:20
load('F:/irdmRound2/stateData.mat')
% Get rid of NaN animals
thisData = cell(1,size(data,2));
c = 1;
for ii = 1:size(data,1)
    if ~any(isnan(data{ii,1})) && ~isnan(data{ii,4})
        thisData(c,:) = data(ii,:);
        c = c+1;
    end
end
% Get base
baseData = thisData(cell2mat(thisData(:,5))==0,:);
uID = unique(baseData(:,2));
for ii = 1:numel(uID)
    baseCatX{ii} = cat(1,baseData{contains(baseData(:,2),uID{ii}),1});
    baseCatY{ii} = cat(1,baseData{contains(baseData(:,2),uID{ii}),4});
end
% Get weights
n = cellfun(@(x) size(x,1),baseCatX);
w = 1./(n/sum(n));
for ii = 1:numel(uID)
    baseCatW{ii} = repmat(w(ii),n(ii),1);
end
% Get sample sizes
mInds = contains(sexData.sex,'M');
fInds = contains(sexData.sex,'F');
mSamp = n(mInds);
fSamp = n(fInds);
%% Build L2O male model
baseCatMaleX = baseCatX(mInds);
baseCatMaleY = baseCatY(mInds);
baseCatMaleYlog = cellfun(@(x) x > median(cat(1,baseCatMaleY{:})),...
    baseCatMaleY,'uniformoutput',0);
% Find pairs of males that have both 0s and 1s
cmbsM = nchoosek(1:numel(baseCatMaleYlog),2);
for ii = 1:size(cmbsM,1)
    chk(ii) = numel(unique(cat(1,baseCatMaleYlog{cmbsM(ii,:)}))) == 2;
end
cmbsM(~chk,:) = [];
% Go through all combinations
for ii = 1:size(cmbsM,1)
    disp(ii)
    testMaleX{ii} = cat(1,baseCatMaleX{cmbsM(ii,:)});
    testMaleY{ii} = cat(1,baseCatMaleYlog{cmbsM(ii,:)});
    others = ~ismember(1:numel(baseCatMaleX),cmbsM(ii,:));
    trainX = cat(1,baseCatMaleX{others});
    trainY = cat(1,baseCatMaleYlog{others});
    % Weight contributions of each rat
    n = cellfun(@(x) size(x,1),baseCatMaleX(others));
    w = 1./(n/sum(n));
    trainW = [];
    for jj = 1:numel(w)
        trainW = [trainW;repmat(w(jj),n(jj),1)];
    end
    % Build model
    mdlM{ii} = fitglm(trainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdlM{ii},testMaleX{ii});
    [thisX,thisY,~,aM(ii)] = perfcurve(testMaleY{ii},pred,1);
    xM(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,50));
    yM(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,50));
    % Permuted
    [~,~,~,aPM(ii)] = perfcurve(testMaleY{ii}(randperm(numel(testMaleY{ii}),...
        numel(testMaleY{ii}))),pred,1);
    % Single features
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'weights',trainW,...
            'distribution','binomial');
        pred = predict(mdl,testMaleX{ii}(:,jj));
        coeffM(ii,jj) = table2array(mdl.Coefficients(2,1));
        [~,~,~,aMS(ii,jj)] = perfcurve(testMaleY{ii},pred,1);
    end
    if ii <= 100
        % Weight contributions of each rat
        n = cellfun(@(x) size(x,1),baseCatMaleX(:));
        w = 1./(n/sum(n));
        trainW = [];
        for jj = 1:numel(w)
            trainW = [trainW;repmat(w(jj),n(jj),1)];
        end
        [trainX80,trainY80,testX80,testY80,trainInds] = trainTest(...
            [trainX;testMaleX{1}],[trainY;testMaleY{1}],0.2);
        mdl80 = fitglm(trainX80,trainY80,'weights',trainW(trainInds),...
            'distribution','binomial');
        pred = predict(mdl80,testX80);
        [~,~,~,a80(ii)] = perfcurve(testY80,pred,1);

        [~,~,~,a80p(ii)] = perfcurve(testY80(randperm(numel(testY80),...
            numel(testY80))),pred,1);
    end
end
%% Build L2O female model; also subset down to 12 rats in training set
baseCatFemaleX = baseCatX(fInds);
baseCatFemaleY = baseCatY(fInds);
baseCatFemaleYlog = cellfun(@(x) x > median(cat(1,baseCatFemaleY{:})),...
    baseCatFemaleY,'uniformoutput',0);
% Find pairs of females that have both 0s and 1s
cmbsF = nchoosek(1:numel(baseCatFemaleYlog),2);
for ii = 1:size(cmbsF,1)
    chk(ii) = numel(unique(cat(1,baseCatFemaleYlog{cmbsF(ii,:)}))) == 2;
end
cmbsF(~chk,:) = [];
% Go through all combinations
for ii = 1:size(cmbsF,1)
    disp(ii)
    testFemaleX{ii} = cat(1,baseCatFemaleX{cmbsF(ii,:)});
    testFemaleY{ii} = cat(1,baseCatFemaleYlog{cmbsF(ii,:)});
    others = logicFind(1,~ismember(1:numel(baseCatFemaleX),cmbsF(ii,:)),...
        '==');
    % Only grab 12
    others12 = randperm(numel(others),12);
    trainX = cat(1,baseCatFemaleX{others12});
    trainY = cat(1,baseCatFemaleYlog{others12});
    % Weight contributions of each rat
    n = cellfun(@(x) size(x,1),baseCatFemaleX(others12));
    w = 1./(n/sum(n));
    trainW = [];
    for jj = 1:numel(w)
        trainW = [trainW;repmat(w(jj),n(jj),1)];
    end
    % Build model
    mdlF{ii} = fitglm(trainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdlF{ii},testFemaleX{ii});
    [thisX,thisY,~,aF(ii)] = perfcurve(testFemaleY{ii},pred,1);
    xF(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,50));
    yF(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,50));
    % Permuted
    [~,~,~,aPF(ii)] = perfcurve(testFemaleY{ii}(randperm(numel(testFemaleY{ii}),...
        numel(testFemaleY{ii}))),pred,1);
    % Feature permutation (beta and high gamma)
    feat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
    hgInds = logicFind(1,contains(feat(1:216),'hg'),'==');
    hgTrainX = trainX;
    for h = 1:numel(hgInds)
        hgTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),h);
    end
    mdl = fitglm(hgTrainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX{ii});
    [~,~,~,aFHG(ii)] = perfcurve(testFemaleY{ii},pred,1);

    bInds = logicFind(1,contains(feat(1:216),'b'),'==');
    bTrainX = trainX;
    for b = 1:numel(bInds)
        bTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),b);
    end
    mdl = fitglm(bTrainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX{ii});
    [~,~,~,aFB(ii)] = perfcurve(testFemaleY{ii},pred,1);

    ofclInds = logicFind(1,contains(feat(1:216),'lOFC'),'==');
    ofclTrainX = trainX;
    for o = 1:numel(ofclInds)
        ofclTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(ofclTrainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX{ii});
    [~,~,~,aFOFCl(ii)] = perfcurve(testFemaleY{ii},pred,1);

    ofcrInds = logicFind(1,contains(feat(1:216),'rOFC'),'==');
    ofcrTrainX = trainX;
    for o = 1:numel(ofcrInds)
        ofcrTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(ofcrTrainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX{ii});
    [~,~,~,aFOFCr(ii)] = perfcurve(testFemaleY{ii},pred,1);

    pfclInds = logicFind(1,contains(feat(1:216),'lPFC'),'==');
    pfclTrainX = trainX;
    for p = 1:numel(pfclInds)
        pfclTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(pfclTrainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX{ii});
    [~,~,~,aFPFCl(ii)] = perfcurve(testFemaleY{ii},pred,1);

    pfcrInds = logicFind(1,contains(feat(1:216),'rPFC'),'==');
    pfcrTrainX = trainX;
    for p = 1:numel(pfcrInds)
        pfcrTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(pfcrTrainX,trainY,'weights',trainW,'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX{ii});
    [~,~,~,aFPFCr(ii)] = perfcurve(testFemaleY{ii},pred,1);
    if ii <= 100 % Run 80:20
        % Weight contributions of each rat
        % Grab 14 random rats
        rats = randperm(numel(baseCatFemaleYlog),17);
        n = cellfun(@(x) size(x,1),baseCatFemaleX(rats));
        w = 1./(n/sum(n));
        trainW = [];
        for jj = 1:numel(w)
            trainW = [trainW;repmat(w(jj),n(jj),1)];
        end
        [trainX80,trainY80,testX80,testY80,trainInds] = trainTest(...
            cat(1,baseCatFemaleX{rats}),cat(1,baseCatFemaleYlog{rats}),0.2);
        mdl80 = fitglm(trainX80,trainY80,'weights',trainW(trainInds),...
            'distribution','binomial');
        pred = predict(mdl80,testX80);
        [~,~,~,a80F(ii)] = perfcurve(testY80,pred,1);

        [~,~,~,aP80F(ii)] = perfcurve(testY80(randperm(numel(testY80),...
            numel(testY80))),pred,1);
    end
    % Single features
%     for jj = 1:216
%         mdl = fitglm(trainX(:,jj),trainY,'weights',trainW,...
%             'distribution','binomial');
%         pred = predict(mdl,testFemaleX{ii}(:,jj));
%         coeffF(ii,jj) = table2array(mdl.Coefficients(2,1));
%         [~,~,~,aFS(ii,jj)] = perfcurve(testFemaleY{ii},pred,1);
%     end
end
%% Compare LOO performance to rats held out
figure
hold on
for ii = 1:17
    plot(ii,aF(any(cmbsF==ii,2)),'.k')
    plot(ii,mean(aF(any(cmbsF==ii,2))),'or')
end
%% 80:20
baseCatFemaleX = baseCatX(fInds);
baseCatFemaleY = baseCatY(fInds);
baseCatFemaleYlog = cellfun(@(x) x > median(cat(1,baseCatFemaleY{:})),...
    baseCatFemaleY,'uniformoutput',0);
for ii = 1:100
    disp(ii)
    inds12 = randperm(numel(baseCatFemaleX),12);
    [trainX,trainY,testFemaleX,testFemaleY,inds] = trainTest(cat(1,baseCatFemaleX{inds12}),...
        cat(1,baseCatFemaleYlog{inds12}),0.2);
    % Weight contributions of each rat
    n = cellfun(@(x) size(x,1),baseCatFemaleX(inds12));
    w = 1./(n/sum(n));
    trainW = [];
    for jj = 1:numel(w)
        trainW = [trainW;repmat(w(jj),n(jj),1)];
    end
    % Build model
    mdlF{ii} = fitglm(trainX,trainY,'weights',trainW(inds),'distribution',...
        'binomial');
    pred = predict(mdlF{ii},testFemaleX);
    [thisX,thisY,~,aF(ii)] = perfcurve(testFemaleY,pred,1);
    xF(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,50));
    yF(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,50));
    % Permuted
    [~,~,~,aPF(ii)] = perfcurve(testFemaleY(randperm(numel(testFemaleY),...
        numel(testFemaleY))),pred,1);
    % Feature permutation (beta and high gamma)
    feat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
        'rNAcC'},{'d','t','a','b','lg','hg'});
    hgInds = logicFind(1,contains(feat(1:216),'hg'),'==');
    hgTrainX = trainX;
    for h = 1:numel(hgInds)
        hgTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),h);
    end
    mdl = fitglm(hgTrainX,trainY,'weights',trainW(inds),'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX);
    [~,~,~,aFHG(ii)] = perfcurve(testFemaleY,pred,1);

    bInds = logicFind(1,contains(feat(1:216),'b'),'==');
    bTrainX = trainX;
    for b = 1:numel(bInds)
        bTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),b);
    end
    mdl = fitglm(bTrainX,trainY,'weights',trainW(inds),'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX);
    [~,~,~,aFB(ii)] = perfcurve(testFemaleY,pred,1);

    ofclInds = logicFind(1,contains(feat(1:216),'lOFC'),'==');
    ofclTrainX = trainX;
    for o = 1:numel(ofclInds)
        ofclTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(ofclTrainX,trainY,'weights',trainW(inds),'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX);
    [~,~,~,aFOFCl(ii)] = perfcurve(testFemaleY,pred,1);

    ofcrInds = logicFind(1,contains(feat(1:216),'rOFC'),'==');
    ofcrTrainX = trainX;
    for o = 1:numel(ofcrInds)
        ofcrTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(ofcrTrainX,trainY,'weights',trainW(inds),'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX);
    [~,~,~,aFOFCr(ii)] = perfcurve(testFemaleY,pred,1);

    pfclInds = logicFind(1,contains(feat(1:216),'lPFC'),'==');
    pfclTrainX = trainX;
    for p = 1:numel(pfclInds)
        pfclTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(pfclTrainX,trainY,'weights',trainW(inds),'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX);
    [~,~,~,aFPFCl(ii)] = perfcurve(testFemaleY,pred,1);

    pfcrInds = logicFind(1,contains(feat(1:216),'rPFC'),'==');
    pfcrTrainX = trainX;
    for p = 1:numel(pfcrInds)
        pfcrTrainX(:,h) = trainX(randperm(size(trainX,1),size(trainX,1)),o);
    end
    mdl = fitglm(pfcrTrainX,trainY,'weights',trainW(inds),'distribution',...
        'binomial');
    pred = predict(mdl,testFemaleX);
    [~,~,~,aFPFCr(ii)] = perfcurve(testFemaleY,pred,1);
end
%% Plot M and F models as histograms with pvalues and permuted
figure
hold on
[f,xi,bw] = ksdensity(aPF);
plot(xi,f*bw)
[f,xi,bw] = ksdensity(aF,'BandWidth',0.05);
plot(xi,f*bw)
plot([mean(aFB) mean(aFB)],[0 0.25])
plot([mean(aFHG) mean(aFHG)],[0 0.25])
plot([mean(aFOFCl) mean(aFOFCl)],[0 0.25])
plot([mean(aFPFCl) mean(aFPFCl)],[0 0.25])
plot([mean(aFOFCr) mean(aFOFCr)],[0 0.25])
plot([mean(aFPFCr) mean(aFPFCr)],[0 0.25])
legend({['permuted: ',num2str(round(mean(aPF),2)),'\pm',...
    num2str(round(conf(aPF,0.95),2))],...
    ['real: ',num2str(round(mean(aF),2)),'\pm',...
    num2str(round(conf(aF,0.95),2)),' ',...
    num2str((sum(aPF>mean(aF))+1)/(numel(aPF)+1))],...
    ['B permute: ',num2str(round(mean(aFB),2)),'\pm',...
    num2str(round(conf(aFB,0.95),2)),' ',...
    num2str((sum(aF<mean(aFB))+1)/(numel(aF)+1))],...
    ['HG permute: ',num2str(round(mean(aFHG),2)),'\pm',...
    num2str(round(conf(aFHG,0.95),2)),' ',...
    num2str((sum(aF<mean(aFHG))+1)/(numel(aF)+1))],...
    ['OFC l permute: ',num2str(round(mean(aFOFCl),2)),'\pm',...
    num2str(round(conf(aFOFCl,0.95),2)),' ',...
    num2str((sum(aF<mean(aFOFCl))+1)/(numel(aF)+1))],...
    ['IL l permute: ',num2str(round(mean(aFPFCl),2)),'\pm',...
    num2str(round(conf(aFPFCl,0.95),2)),' ',...
    num2str((sum(aF<mean(aFPFCl))+1)/(numel(aF)+1))],...
    ['OFC r permute: ',num2str(round(mean(aFOFCr),2)),'\pm',...
    num2str(round(conf(aFOFCr,0.95),2)),' ',...
    num2str((sum(aF<mean(aFOFCr))+1)/(numel(aF)+1))],...
    ['IL r permute: ',num2str(round(mean(aFPFCr),2)),'\pm',...
    num2str(round(conf(aFPFCr,0.95),2)),' ',...
    num2str((sum(aF<mean(aFPFCr))+1)/(numel(aF)+1))]})
%% Male
figure
hold on
[f,xi,bw] = ksdensity(aPM,'BandWidth',0.05);
plot(xi,f*bw)
plot([mean(aM) mean(aM)],[0 0.25])
legend({['permuted: ',num2str(round(mean(aPM),2)),'\pm',...
    num2str(round(conf(aPM,0.95),2))],...
    ['real: ',num2str(round(mean(aM),2)),'\pm',...
    num2str(round(conf(aM,0.95),2)),' ',...
    num2str((sum(aPM>mean(aM))+1)/(numel(aPM)+1))]})
%% Male and female
figure
hold on
[f,xi,bw] = ksdensity(aPM,'BandWidth',0.05);
plot(xi,f*bw)
plot([mean(a80) mean(a80)],[0 0.25])
[f,xi,bw] = ksdensity(aPF);
plot(xi,f*bw)
plot([mean(a80F) mean(a80F)],[0 0.25])
legend({['permuted M: ',num2str(round(mean(aPM),2)),'\pm',...
    num2str(round(conf(aPM,0.95),2))],...
    ['M: ',num2str(round(mean(a80),2)),'\pm',...
    num2str(round(conf(a80,0.95),2)),' ',...
    num2str((sum(aPM>mean(a80))+1)/(numel(aPM)+1))],...
    ['permuted F: ',num2str(round(mean(aPF),2)),'\pm',...
    num2str(round(conf(aPF,0.95),2))],...
    ['F: ',num2str(round(mean(a80F),2)),'\pm',...
    num2str(round(conf(a80F,0.95),2)),' ',...
    num2str((sum(aPF>mean(a80F))+1)/(numel(aPF)+1))]})
%% Cross apply models - way more iterations by testing each test set on 
% each model...
for ii = 1:numel(mdlF)
    for jj = 1:numel(testMaleX)
        pred = predict(mdlF{ii},testMaleX{jj});
        [thisX,thisY,~,aFM(ii,jj)] = perfcurve(testMaleY{jj},pred,1);
        xFM(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
            linspace(0,1,50));
        yFM(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
            linspace(0,1,50));
    end
end
for ii = 1:numel(mdlM)
    for jj = 1:numel(testFemaleX)
        pred = predict(mdlM{ii},testFemaleX{jj});
        [thisX,thisY,~,aMF(ii,jj)] = perfcurve(testFemaleY{jj},pred,1);
        xMF(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
            linspace(0,1,50));
        yMF(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
            linspace(0,1,50));
    end
end
% save('F:\irdmRound2\ddtStateSex.mat','xM','yM','aM','coeffM','aMS',...
%     'mdlM','xF','yF','aF','aFM','coeffF','aFS','xFM','yFM','xMF','yMF')
%% Plot M and F models
figure
subplot(1,2,1)
hold on
plot(mean(xM,1),mean(yM))
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')

plot(mean(xF,1),mean(yF))
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
legend({['male: ',num2str(round(mean(aM),2)),'\pm',...
    num2str(round(conf(aM,0.95),2))],['female: ',...
    num2str(round(mean(aF),2)),'\pm',num2str(round(conf(aF,0.95),2))]})

subplot(1,2,2)
hold on
plot(squeeze(mean(xFM,[1,2])),squeeze(mean(yFM,[1,2])))
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')

plot(squeeze(mean(xMF,[1,2])),squeeze(mean(yMF,[1,2])))
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
legend({['f->m: ',num2str(round(mean(reshape(aFM,1,numel(aFM))),2)),'\pm',...
    num2str(round(conf(reshape(aFM,1,numel(aFM)),0.95),2))],['m->f: ',...
    num2str(round(mean(reshape(aMF,1,numel(aMF))),2)),'\pm',...
    num2str(round(conf(reshape(aMF,1,numel(aMF)),0.95),2))]})

figure
subplot(1,2,1)
histogram(cat(1,baseCatMaleY{:}),'normalization','probability')
title('male AUC')
hold on
plot([median(cat(1,baseCatMaleY{:})) median(cat(1,baseCatMaleY{:}))],...
    [0 1],'--k')
subplot(1,2,2)
histogram(cat(1,baseCatFemaleY{:}),'normalization','probability')
title('female AUC')
hold on
plot([median(cat(1,baseCatFemaleY{:})) median(cat(1,baseCatFemaleY{:}))],...
    [0 1],'--k')
%% Compare single features - male vs. female
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
[aFSsort,sortInds] = sort(mean(aFS,1),'descend');
fFeat = feat(sortInds)';
f = mean(aFS,1).*sign(mean(coeffF,1));

[aMSsort,sortInds] = sort(mean(aMS,1),'descend');
mFeat = feat(sortInds)';
m = mean(aMS,1).*sign(mean(coeffM,1));
figure
plot(f,m,'ok')
%% Load lasso results (male)
cd('F:\irdmRound2\stateDDT\')
for ii = 1:87
    load(['ddtStateM',num2str(ii),'.mat'])
    aM(ii) = accM{1}.acc;
    xM(ii,:) = interp1(linspace(0,1,numel(accM{1}.x)),accM{1}.x,...
            linspace(0,1,50));
    yM(ii,:) = interp1(linspace(0,1,numel(accM{1}.y)),accM{1}.y,...
            linspace(0,1,50));
end
for ii = 1:134
    load(['ddtStateF',num2str(ii),'.mat'])
    aF(ii) = accF{1}.acc;
    xF(ii,:) = interp1(linspace(0,1,numel(accF{1}.x)),accF{1}.x,...
            linspace(0,1,50));
    yF(ii,:) = interp1(linspace(0,1,numel(accF{1}.y)),accF{1}.y,...
            linspace(0,1,50));
    
end

figure
hold on
plot(mean(xM),mean(yM))
plot(mean(xF),mean(yF))
%% Build combined male & female model
% still total of 12 in training and two in test, so 6 male and 6 female in
% train and one of each in test

% Redetermine median based on whole population
baseCatFemaleYlogCmb = cellfun(@(x) x > median(cat(1,baseCatFemaleY{:},baseCatMaleY{:})),...
    baseCatFemaleY,'uniformoutput',0);
baseCatMaleYlogCmb = cellfun(@(x) x > median(cat(1,baseCatFemaleY{:},baseCatMaleY{:})),...
    baseCatMaleY,'uniformoutput',0);
% Find combinations of males and females that have 1s and 0s for test
c = 1;
for ii = 1:numel(baseCatMaleYlog)
    for jj = 1:numel(baseCatFemaleYlog)
        if numel(unique(cat(1,baseCatMaleYlogCmb{ii},...
                baseCatFemaleYlogCmb{jj}))) == 2
            cmbs(c,:) = [ii,jj];
            c = c+1;
        end
    end
end
predM = []; testYM = []; predF = []; testYF = [];
for ii = 1:size(cmbs,1)
    testY =  cat(1,baseCatMaleYlogCmb{cmbs(ii,1)},...
        baseCatFemaleYlogCmb{cmbs(ii,2)});
    testX = cat(1,baseCatMaleX{cmbs(ii,1)},...
        baseCatFemaleX{cmbs(ii,2)});
    % Find 6 male and females of the others for training
    othersF = logicFind(1,~ismember(1:numel(baseCatFemaleX),cmbs(ii,2)),...
        '==');
    othersM = logicFind(1,~ismember(1:numel(baseCatMaleX),cmbs(ii,1)),...
        '==');
    % Only grab 6
    others6F = randperm(numel(othersF),6);
    others6M = randperm(numel(othersM),6);
    trainX = cat(1,baseCatMaleX{others6M},baseCatFemaleX{others6F});
    trainY = cat(1,baseCatMaleYlogCmb{others6M},...
        baseCatFemaleYlogCmb{others6F});
    % Weight contributions of each rat
    n = cat(2,cellfun(@(x) size(x,1),baseCatMaleX(others6M)),...
        cellfun(@(x) size(x,1),baseCatFemaleX(others6F)));
    w = 1./(n/sum(n));
    trainW = [];
    for jj = 1:numel(w)
        trainW = [trainW;repmat(w(jj),n(jj),1)];
    end
    % Build model
    mdl = fitglm(trainX,trainY,'weights',trainW,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,aCmb(ii)] = perfcurve(testY,pred,1);
    xCmb(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
    yCmb(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));

    [~,~,~,aCmbP(ii)] = perfcurve(testY(randperm(numel(testY))),pred,1);
    for jj = 1:216
        smdl = fitglm(trainX(:,jj),trainY,'weights',trainW,...
            'distribution','binomial');
        pred = predict(smdl,testX(:,jj));
        [~,~,~,aCmbS(ii,jj)] = perfcurve(testY,pred,1);
        cmbCoeff(ii,jj) = table2array(smdl.Coefficients(2,1));
    end
    % Test on male - will need extra animal, which makes interpretation
    % harder
    leftover = logicFind(0,ismember(othersM,others6M),'==');
    % Pick random leftover
    testYM = [];
    while numel(unique(testYM)) ~= 2
        ind = randi(numel(leftover),1);
        testYM = cat(1,baseCatMaleYlogCmb{cmbs(ii,1)},...
            baseCatMaleYlogCmb{leftover(ind)});
        testXM = cat(1,baseCatMaleX{cmbs(ii,1)},...
            baseCatMaleX{leftover(ind)});
    end    
    pred = predict(mdl,testXM);
    [thisX,thisY,~,aCmbM(ii)] = perfcurve(...
        testYM,pred,1);
    xCmbM(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
    yCmbM(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));
%     % Store prediction on male
%     predM = [predM;predict(mdl,baseCatMaleX{cmbs(ii,1)})];
%     testYM = [testYM;baseCatMaleYlog{cmbs(ii,1)}];

    % Test on female - will need extra animal
    leftover = logicFind(0,ismember(othersF,others6F),'==');
    % Pick random leftover
    testYF = [];
    while numel(unique(testYF)) ~= 2
        ind = randi(numel(leftover),1);
        testYF = cat(1,baseCatFemaleYlogCmb{cmbs(ii,2)},...
            baseCatFemaleYlogCmb{leftover(ind)});
        testXF = cat(1,baseCatFemaleX{cmbs(ii,2)},...
            baseCatFemaleX{leftover(ind)});
    end
    pred = predict(mdl,testXF);
    [thisX,thisY,~,aCmbF(ii)] = perfcurve(testYF,pred,1);
    xCmbF(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
    yCmbF(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));
%     % Store prediction on female
%     predF = [predF;predict(mdl,baseCatFemaleX{cmbs(ii,2)})];
%     testYF = [testYF;baseCatFemaleYlog{cmbs(ii,2)}];
    % Test model on combined 2 male 2 female
    pred = predict(mdl,[testXM;testXF]);
    [thisX,thisY,~,aCmbMF(ii)] = perfcurve([testYM;testYF],pred,1);
    xCmbMF(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
    yCmbMF(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));
    disp(ii)
end
% [xm,ym,~,am] = perfcurve(testYM,predM,1);
% [xf,yf,~,af] = perfcurve(testYF,predF,1);

save('F:\irdmRound2\ddtStateSexCmb.mat','aCmb','xCmb','yCmb','aCmbF',...
    'xCmbF','yCmbF','aCmbM','xCmbM','yCmbM','aCmbMF','xCmbMF','yCmbMF',...
    'aCmbP')
%% Plot
figure
hold on
plot(mean(xCmb,1),mean(yCmb))
plot(mean(xCmbM,1),mean(yCmbM,1))
plot(mean(xCmbF,1),mean(yCmbF,1))
legend({['gen->gen: ',num2str(round(mean(aCmb),2)),'\pm',...
    num2str(round(conf(aCmb,0.95),2))],['gen->m: ',num2str(round(mean(aCmbM),2)),'\pm',...
    num2str(round(conf(aCmbM,0.95),2))],['gen->f: ',...
    num2str(round(mean(aCmbF),2)),'\pm',num2str(round(conf(aCmbF,0.95)...
    ,2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
%% histograms and pvalues
figure
hold on
[f,xi,bw] = ksdensity(aCmbP);
plot(xi,f*bw)
plot([mean(aCmb) mean(aCmb)],[0 0.25])
plot([mean(aCmbM) mean(aCmbM)],[0 0.25])
plot([mean(aCmbF) mean(aCmbF)],[0 0.25])
legend({['permuted: ',num2str(round(mean(aCmbP),2)),'\pm',...
    num2str(round(conf(aCmbP,0.95),2))],...
    ['gen->gen: ',num2str(round(mean(aCmb),2)),'\pm',...
    num2str(round(conf(aCmb,0.95),2)),' ',...
    num2str((sum(aCmbP>mean(aCmb))+1)/(numel(aCmbP)+1))],['gen->M: ',...
    num2str(round(mean(aCmbM),2)),'\pm',...
    num2str(round(conf(aCmbM,0.95),2)),' ',...
    num2str((sum(aCmbP>mean(aCmbM))+1)/(numel(aCmbP)+1))],['gen->F: ',...
    num2str(round(mean(aCmbF),2)),'\pm',...
    num2str(round(conf(aCmbF,0.95),2)),' ',...
    num2str((sum(aCmbP>mean(aCmbF))+1)/(numel(aCmbP)+1))]})
set(gca,'XLim',[0 1],'ylim',[0 0.2])
%% Compare single features - female vs. cmb
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
% load('F:/irdmRound2/ddtStateSex.mat','aFS','coeffF')
[~,sortInds] = sort(mean(aCmbS,1)','descend');
cmbFeat = feat(sortInds)';
cmb = mean(aCmbS,1).*sign(mean(cmbCoeff,1));
aCmbSort = cmb(sortInds)';
% load('F:/irdmRound2/ddtStateSexCmb.mat')
[~,sortInds] = sort(mean(aFS,1),'descend');
fFeat = feat(sortInds)';
f = mean(aFS,1).*sign(mean(coeffF,1));
aFSort = f(sortInds)';
c = distinguishable_colors(6);
figure
limx = [-0.8 -0.4;0.4 0.8;-0.8 -0.4;0.4 0.8];
limy = [0.4 0.8;0.4 0.8;-0.8 -0.4;-0.8 -0.4];
for ii = 1:4
    subplot(2,2,ii)
    hold on
    scatter(f(1:48),cmb(1:48),[],repmat(c,8,1),'o');
    scatter(f(49:end),cmb(49:end),[],repmat(c,28,1),'s');
    set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
    xlim(limx(ii,1:2))
    ylim(limy(ii,1:2))
end
%% Build combined male & female continuous model
% still total of 12 in training and two in test, so 6 male and 6 female in
% train and one of each in test
load('F:/irdmRound2/stateData.mat')
% Get rid of NaN animals
thisData = cell(1,size(data,2));
c = 1;
for ii = 1:size(data,1)
    if ~any(isnan(data{ii,1})) && ~isnan(data{ii,4})
        thisData(c,:) = data(ii,:);
        c = c+1;
    end
end
% Get base
baseData = thisData(cell2mat(thisData(:,5))==0,:);
uID = unique(baseData(:,2));
for ii = 1:numel(uID)
    baseCatX{ii} = cat(1,baseData{contains(baseData(:,2),uID{ii}),1});
    baseCatY{ii} = cat(1,baseData{contains(baseData(:,2),uID{ii}),4});
end
% Get weights
n = cellfun(@(x) size(x,1),baseCatX);
w = 1./(n/sum(n));
for ii = 1:numel(uID)
    baseCatW{ii} = repmat(w(ii),n(ii),1);
end
% Get male and female data
mInds = contains(sexData.sex,'M');
fInds = contains(sexData.sex,'F');
mSamp = n(mInds);
fSamp = n(fInds);
baseCatFemaleX = baseCatX(fInds);
baseCatFemaleY = baseCatY(fInds);
baseCatMaleX = baseCatX(mInds);
baseCatMaleY = baseCatY(mInds);
% Find combinations of males and females
c = 1;
for ii = 1:numel(baseCatMaleY)
    for jj = 1:numel(baseCatFemaleY)
        cmbs(c,:) = [ii,jj];
        c = c+1;
    end
end
% Go through all combinations
predM = []; testYM = []; predF = []; testYF = [];
for ii = 1:size(cmbs,1)
%     testY =  cat(1,baseCatMaleY{cmbs(ii,1)},...
%         baseCatFemaleY{cmbs(ii,2)});
%     testX = cat(1,baseCatMaleX{cmbs(ii,1)},...
%         baseCatFemaleX{cmbs(ii,2)});
%     % Find 6 male and females of the others for training
%     othersF = logicFind(1,~ismember(1:numel(baseCatFemaleX),cmbs(ii,2)),...
%         '==');
%     othersM = logicFind(1,~ismember(1:numel(baseCatMaleX),cmbs(ii,1)),...
%         '==');
%     % Only grab 6
%     others6F = randperm(numel(othersF),6);
%     others6M = randperm(numel(othersM),6);
%     trainX = cat(1,baseCatMaleX{others6M},baseCatFemaleX{others6F});
%     trainY = cat(1,baseCatMaleY{others6M},...
%         baseCatFemaleY{others6F});
%     % Weight contributions of each rat
%     n = cat(2,cellfun(@(x) size(x,1),baseCatMaleX(others6M)),...
%         cellfun(@(x) size(x,1),baseCatFemaleX(others6F)));
%     w = 1./(n/sum(n));
%     trainW = [];
%     for jj = 1:numel(w)
%         trainW = [trainW;repmat(w(jj),n(jj),1)];
%     end
%     % Build model
%     mdl = fitglm(trainX,trainY,'weights',trainW,'distribution','poisson');
%     pred = predict(mdl,testX);
%     maeCmb(ii) = mean(abs(pred-testY));
%     % Permuted test
%     pY = testY(randperm(numel(testY)));
%     maeCmbPerm(ii) = mean(abs(pred-pY));
%     for jj = 1:216
%         smdl = fitglm(trainX(:,jj),trainY,'weights',trainW,...
%             'distribution','normal');
%         pred = predict(smdl,testX(:,jj));
%         maeCmbS(ii,jj) = mean(abs(pred-testY));
%         cmbCoeff(ii,jj) = table2array(smdl.Coefficients(2,1));
%     end
%     % Test on male - will need extra animal, which makes interpretation
%     % harder
%     leftover = logicFind(0,ismember(othersM,others6M),'==');
%     % Pick random leftover
%     ind = randi(numel(leftover),1);
%     testYM = cat(1,baseCatMaleY{cmbs(ii,1)},...
%         baseCatMaleY{leftover(ind)});
%     testXM = cat(1,baseCatMaleX{cmbs(ii,1)},...
%         baseCatMaleX{leftover(ind)});
%     pred = predict(mdl,testXM);
%     maeCmbM(ii) = mean(abs(pred-testYM));
%     % Test on female - will need extra animal
%     leftover = logicFind(0,ismember(othersF,others6F),'==');
%     % Pick random leftover
%     ind = randi(numel(leftover),1);
%     testYF = cat(1,baseCatFemaleY{cmbs(ii,2)},...
%         baseCatFemaleY{leftover(ind)});
%     testXF = cat(1,baseCatFemaleX{cmbs(ii,2)},...
%         baseCatFemaleX{leftover(ind)});
%     pred = predict(mdl,testXF);
%     maeCmbF(ii) = mean(abs(pred-testYF));
%     % Test model on combined 2 male 2 female
%     pred = predict(mdl,[testXM;testXF]);
%     maeCmbMF(ii) = mean(abs(pred-[testYM;testYF]));
%     disp(ii)
    % 80:20 model using 7 of each sex
    male7 = randperm(numel(baseCatMaleX),7);
    female7 = randperm(numel(baseCatFemaleX),7);
    fullX = [cat(1,baseCatMaleX{male7});cat(1,baseCatFemaleX{female7})];
    fullY = [cat(1,baseCatMaleY{male7});cat(1,baseCatFemaleY{female7})];
    % Weight contributions of each rat
    n = cat(2,cellfun(@(x) size(x,1),baseCatMaleX(male7)),...
        cellfun(@(x) size(x,1),baseCatFemaleX(female7)));
    w = 1./(n/sum(n));
    trainW = [];
    for jj = 1:numel(w)
        trainW = [trainW;repmat(w(jj),n(jj),1)];
    end
    [~,~,~,~,trainInd,testInd] = trainTest([1:numel(fullY)]',1:numel(fullY),0.2);
    testY = fullY(testInd);
    mdl = fitglm(fullX(trainInd,:),fullY(trainInd),...
        'weights',trainW(trainInd),'distribution','poisson');
    pred = predict(mdl,fullX(testInd,:));
    maeCmb8020(ii) = mean(abs(pred-testY));
    % Permuted test
    pY = testY(randperm(numel(testY)));
    maeCmbPerm(ii) = mean(abs(pred-pY));
end
% save('F:\irdmRound2\ddtStateSexCmbCont.mat','maeCmb','maeCmbF','maeCmbM'...
%     ,'maeCmbMF','maeCmbS')
%%
figure
h(1) = histogram(maeCmb,'binwidth',5,'normalization','probability');
hold on
h(2) = plot([mean(maeCmb) mean(maeCmb)],[0 0.4],'b');
h(3) = histogram(maeCmbPerm,'binwidth',5,'normalization','probability');
h(4) = plot([mean(maeCmbPerm) mean(maeCmbPerm)],[0 0.4],'r');
box off
xlabel('MAE')
ylabel('proportion of models')
legend(h([2 4]),{['real: ',num2str(round(mean(maeCmb),2)),'\pm',...
    num2str(round(conf(maeCmb,0.95),2))],...
    ['perm: ',num2str(round(mean(maeCmbPerm),2)),'\pm',...
    num2str(round(conf(maeCmbPerm,0.95),2))]})
title(['generalized; train 6M&6F; test 1M&1F; cmbs = ',...
    num2str(numel(maeCmb))])
figure
hold on
histogram(maeCmbM,'binwidth',5,'normalization','probability')
h(1) = plot([mean(maeCmbM) mean(maeCmbM)],[0 0.4],'b');
histogram(maeCmbF,'binwidth',5,'normalization','probability')
h(2) = plot([mean(maeCmbF) mean(maeCmbF)],[0 0.4],'r');
xlabel('MAE')
ylabel('proportion of models')
legend(h,{['male: ',num2str(round(mean(maeCmbM),2)),'\pm',...
    num2str(round(conf(maeCmbM,0.95),2))],...
    ['female: ',num2str(round(mean(maeCmbF),2)),'\pm',...
    num2str(round(conf(maeCmbF,0.95),2))]})
title(['gen->sex; train 6M&6F; test 2M|2F; cmbs = ',num2str(numel(maeCmb))])
%% Build L2O male model continuous
baseCatMaleX = baseCatX(mInds);
baseCatMaleY = baseCatY(mInds);
% Find pairs of males
cmbsM = nchoosek(1:numel(baseCatMaleY),2);
% Go through all combinations
for ii = 1:size(cmbsM,1)
    disp(ii)
    testMaleX{ii} = cat(1,baseCatMaleX{cmbsM(ii,:)});
    testMaleY{ii} = cat(1,baseCatMaleY{cmbsM(ii,:)});
    others = ~ismember(1:numel(baseCatMaleX),cmbsM(ii,:));
    trainX = cat(1,baseCatMaleX{others});
    trainY = cat(1,baseCatMaleY{others});
    % Weight contributions of each rat
    n = cellfun(@(x) size(x,1),baseCatMaleX(others));
    w = 1./(n/sum(n));
    trainW = [];
    for jj = 1:numel(w)
        trainW = [trainW;repmat(w(jj),n(jj),1)];
    end
    % Build model
    mdlM{ii} = fitglm(trainX,trainY,'weights',trainW,'distribution',...
        'poisson');
    pred = predict(mdlM{ii},testMaleX{ii});
    maeM(ii) = mean(abs(pred-testMaleY{ii}));
    % Single features
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'weights',trainW,...
            'distribution','poisson');
        pred = predict(mdl,testMaleX{ii}(:,jj));
        coeffM(ii,jj) = table2array(mdl.Coefficients(2,1));
        maeMS(ii,jj) = mean(abs(pred-testMaleY{ii}));
    end
    % Permuted
    pY = testMaleY{ii}(randperm(numel(testMaleY{ii})));
    maePermM(ii) = mean(abs(pred-pY)); 
end
% Build L2O female model continuous; also subset down to 12 rats in training set
baseCatFemaleX = baseCatX(fInds);
baseCatFemaleY = baseCatY(fInds);
% Find pairs of females
cmbsF = nchoosek(1:numel(baseCatFemaleY),2);
% Go through all combinations
for ii = 1:size(cmbsF,1)
    disp(ii)
    testFemaleX{ii} = cat(1,baseCatFemaleX{cmbsF(ii,:)});
    testFemaleY{ii} = cat(1,baseCatFemaleY{cmbsF(ii,:)});
    others = logicFind(1,~ismember(1:numel(baseCatFemaleX),cmbsF(ii,:)),...
        '==');
    % Only grab 12
    others12 = randperm(numel(others),12);
    trainX = cat(1,baseCatFemaleX{others12});
    trainY = cat(1,baseCatFemaleY{others12});
    % Weight contributions of each rat
    n = cellfun(@(x) size(x,1),baseCatFemaleX(others12));
    w = 1./(n/sum(n));
    trainW = [];
    for jj = 1:numel(w)
        trainW = [trainW;repmat(w(jj),n(jj),1)];
    end
    % Build model
    mdlF{ii} = fitglm(trainX,trainY,'weights',trainW,'distribution',...
        'poisson');
    pred = predict(mdlF{ii},testFemaleX{ii});
    maeF(ii) = mean(abs(pred-testFemaleY{ii}));
    % Single features
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'weights',trainW,...
            'distribution','poisson');
        pred = predict(mdl,testFemaleX{ii}(:,jj));
        coeffF(ii,jj) = table2array(mdl.Coefficients(2,1));
        maeFS(ii,jj) = mean(abs(pred-testFemaleY{ii}));
    end
    % Permuted
    pY = testFemaleY{ii}(randperm(numel(testFemaleY{ii})));
    maePermF(ii) = mean(abs(pred-pY)); 
end
% save('F:/irdmRound2/ddtStateSexCont.mat','maeM','maeMS','coeffM',...
%     'maeF','maeFS','coeffF')
%%
figure
hold on
histogram(maeM,'binwidth',5,'normalization','probability')
h(1) = plot([mean(maeM) mean(maeM)],[0 1],'b');
histogram(maePermM,'binwidth',5,'normalization','probability')
h(2) = plot([mean(maePermM) mean(maePermM)],[0 1],'r');
xlabel('MAE')
ylabel('proportion of models')
legend(h,{['male: ',num2str(round(mean(maeM),2)),'\pm',...
    num2str(round(conf(maeM,0.95),2))],...
    ['perm: ',num2str(round(mean(maePermM),2)),'\pm',...
    num2str(round(conf(maePermM,0.95),2)),]})
title(['train 12:test 2; male; cmbs = ',num2str(numel(maeM))])
clear h
figure
hold on
histogram(maeF,'binwidth',5,'normalization','probability')
h(1) = plot([mean(maeF) mean(maeF)],[0 1],'b');
histogram(maePermF,'binwidth',5,'normalization','probability')
h(2) = plot([mean(maePermF) mean(maePermF)],[0 1],'r');
xlabel('MAE')
ylabel('proportion of models')
legend(h,{['female: ',num2str(round(mean(maeF),2)),...
    '\pm',num2str(round(conf(maeF,0.95),2))],...
    ['perm: ',num2str(round(mean(maePermF),2)),...
    '\pm',num2str(round(conf(maePermF,0.95),2))]})
title(['train 12:test 2; female; cmbs = ',num2str(numel(maeF))])
%%
% Build 80:20 model
allX = cat(1,baseCatX{:});
allY = cat(1,baseCatY{:});
allW = cat(1,baseCatW{:});
[trainXind,trainYind,testXind,testYind] = trainTest([1:numel(allY)]',allY,0.2);
% mdl = fitglm(allX(trainXind,:),allY(trainXind),'weights',allW(trainXind));
mdl = fitglm(allX(trainXind,:),allY(trainXind));
p = predict(mdl,allX(testXind,:));
mseAll = sum((allY(testXind)-p).^2)/numel(allY(testXind));
%% Median split model
med = median(cell2mat(baseData(:,4)));
allX = cat(1,baseCatX{:});
allY = cat(1,baseCatY{:});
allW = cat(1,baseCatW{:});
[trainXind,trainYind,testXind,testYind] = trainTest([1:numel(allY)]',allY,0.2);
% mdl = fitglm(allX(trainXind,:),allY(trainXind)>=med,'weights',allW(trainXind)...
%     ,'distribution','binomial');
mdl = fitglm(allX(trainXind,:),allY(trainXind)>=med,'distribution',...
    'binomial');
p = predict(mdl,allX(testXind,:));
[~,~,~,a] = perfcurve(allY(testXind)>=med,p,1);