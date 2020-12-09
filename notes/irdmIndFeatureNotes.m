load('irdmDelayImmediateDataMinus1.mat')
% Run individual animal models
allData = squeeze(struct2cell(data));
for n = 1:100
    delay = []; immediate = [];
    disp(n)
    rng(n)
    % Find minimal number of delays and immediates each session contributes
    for jj = 1:size(animalInds,1)
        % Find indices of delays and immediates
        delayInd = logicFind(1,cell2mat(allData(7,animalInds(jj,1):...
            animalInds(jj,2))),'==');
        immediateInd = logicFind(0,cell2mat(allData(7,animalInds(jj,1):...
            animalInds(jj,2))),'==');
        % Count number of delays and immediates
        dN(jj) = numel(delayInd);
        iN(jj) = numel(immediateInd);
        % Find min
        thisMin = min([dN(jj),iN(jj)]);
        if jj == 1
            start = 0;
        else
            start = animalInds(jj-1,2);
        end
        % Pull that number of both trials
        delay = cat(2,delay,allData(:,start+delayInd(...
            randperm(dN(jj),thisMin))));
        immediate = cat(2,immediate,allData(:,start+immediateInd(...
            randperm(iN(jj),thisMin))));
    end
    u = unique(delay(1,:));
    powN = numel(delay{4}(:,:,1,1));
    cohN = numel(delay{6}(:,:,1,1));
    % For combined model with IRDM 11,15,16,22
    theseInds = zeros(1,numel(delay,2));
    for ii = [1,3,4,7]
       theseInds = logical(theseInds+strcmp(u(ii),delay(1,:))); 
    end
    ii = 1;
%     for ii = 1:numel(u)
%         theseInds = strcmp(u(ii),delay(1,:));
        % Delay
        delayRelPower = cellfun(@(x) reshape(squeeze(mean(x,4,...
            'omitnan')),1,powN),delay(4,theseInds),'UniformOutput',0);
        catDelayPower = cat(1,delayRelPower{:});
        delayNormCoh = cellfun(@(x) reshape(squeeze(mean(x,4,...
            'omitnan')),1,cohN),delay(6,theseInds),'UniformOutput',0);
        catDelayCoh = cat(1,delayNormCoh{:});
        % Immediate
        immediateRelPower = cellfun(@(x) reshape(squeeze(mean(x,4,...
            'omitnan')),1,powN),immediate(4,theseInds),'UniformOutput',0);
        catImmediatePower = cat(1,immediateRelPower{:});
        immediateNormCoh = cellfun(@(x) reshape(squeeze(mean(x,4,...
            'omitnan')),1,cohN),immediate(6,theseInds),'UniformOutput',0);
        catImmediateCoh = cat(1,immediateNormCoh{:});
        % Combine
        indDataX{ii} = [catDelayPower,catDelayCoh;catImmediatePower,...
            catImmediateCoh];
        indDataY{ii} = [ones(numel(delayRelPower),1);...
            zeros(numel(immediateRelPower),1)];
        % Split into train and test
        x = indDataX{ii};
        y = indDataY{ii};
        trainInds = randperm(numel(y),ceil(numel(y)*0.8));
        trainX = x(trainInds,:);
        trainY = y(trainInds,:);
        testInds = 1:numel(y);
        testInds = ~ismember(testInds,trainInds);
        testX = x(testInds,:);
        testY = y(testInds);
        % Build models
        for jj = 1:size(trainX,2)
            mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial',...
                'binomialsize',size(trainX,1));
            pred = predict(mdl,testX(:,jj));
            [~,~,~,a(n,ii,jj)] = perfcurve(testY,pred,1);
            % Also grab sign of model; 1 == positive slope; delay is higher
            s(n,ii,jj) = mdl.Coefficients.Estimate(2)/...
                abs(mdl.Coefficients.Estimate(2));
        end
        %         save(['/ihome/ldwiel/data/delayImmediateSingleFeature',u{ii},'_',num2str(n+20*(k-1)),'.mat'],...
        %             'a')
%     end
end

% Full features
fullFeat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
% ch. 1 removed features
minus1Feat = names({'rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
mA = squeeze(mean(a.*s,1));
for ii = 1:7
    [smA(:,ii),ind(ii,:)] = sort(abs(mA(ii,:)),'descend');
    smA(:,ii) = smA(:,ii).*(mA(ii,ind(ii,:))./abs(mA(ii,ind(ii,:))))';
    feat(:,ii) = minus1Feat(ind(ii,:));
end
save('delayImmediateSingleFeatureIndModels.mat','a','s','mA','smA','feat')
%%
eucD = pdist(mA,'euclidean');
clustTreeEuc = linkage(eucD,'average');
figure
dendrogram(clustTreeEuc)
%%
[~,diff25] = sort(abs(mA(2,:) - mA(5,:)),'ascend');
%% Median split
load('G:\GreenLab\data\irdm\ddtDataRaw.mat')
thisY = y;
thisY(y>=median(y)) = 1;
thisY(y<median(y)) = 0;
y = thisY;

% y = [ones(23,1);zeros(12,1);ones(13,1);zeros(32,1);];
for ii = 1:100
    disp(ii)
    rng(ii);
    trainInds = randperm(size(x,1),round(size(x,1)*0.8));
    testInds = ~ismember(1:80,trainInds);
    trainX =  x(trainInds,:);
    trainY = y(trainInds);
    testX = x(testInds,:);
    testY = y(testInds);
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial',...
            'binomialsize',size(trainX,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,a(ii,jj)] = perfcurve(testY,pred,1);
        % Also grab sign of model; 1 == positive slope; 'highs' are higher
        s(ii,jj) = mdl.Coefficients.Estimate(2)/...
            abs(mdl.Coefficients.Estimate(2));
    end
end
%%
% Full features
fullFeat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
% ch. 1 removed features
minus1Feat = names({'rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
mA = squeeze(mean(a.*s,1));
[smA,ind] = sort(abs(mA),'descend');
smA = smA'.*(mA(ind)./abs(mA(ind)))';
feat = fullFeat(ind)';
%% top 4 animal model
for ii = 1:100
    load(['G:\GreenLab\data\irdm\delayImmediateGenTop4\'...
        'delayImmediateGen',num2str(ii),'.mat'])
    [x(ii,:),y(ii,:),~,a(ii)] = perfcurve(hist.cfg.naive.testY,...
        accArray{1}.pred,1,'TVals',0:0.001:1,'UseNearest','off');
end
%%
indA = [];
for ii = 1:100
    load(['G:\GreenLab\data\irdm\delayImmediateGenTop4Ind\'...
        'delayImmediateGen',num2str(ii),'.mat'])
    [x(ii,:),y(ii,:),~,a(ii)] = perfcurve(hist.cfg.naive.testY,...
        accArray{1}.pred,1,'TVals',0:0.001:1,'UseNearest','off');
    indA(:,:,ii) = aInd;
end