load('D:/emily/rip/emilyNVHLdataTrlRaw.mat')
%% Average trials together to minimize variance in data
% avgData = cell(size(data{1},1),1);
% for ii = 1:size(data{1},1)
%     x = 1;
%     while x+29 <= size(data{1}{ii},1) 
%         avgData{ii} = [avgData{ii};mean(data{1}{ii}(x:x+29,:),1)];
%         x = x+29;
%     end
% end
% oldData = data;
% data{1} = avgData;
%% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),data{1,1}(:,1),'uniformoutput',0));
nFiles = size(data{1},1);
[animal,day] = deal(cell(nFiles,1));
for ii = 1:nFiles
    parts = strsplit(files{1}{ii},'_');
    animal{ii,1} = repmat({parts{1}(4:end)},n(ii,1),1);
    day{ii,1} = repmat({parts{2}},n(ii,1),1); %#ok
end
%% Get indices of groups and combine data
SSBind = [1:4,14:17,39:41];
SSB = cat(1,data{1}{SSBind});
SSBid = cat(1,animal{SSBind});

SSIind = [5:6,18:19,42:43];
SSI = cat(1,data{1}{SSIind});
SSIid = cat(1,animal{SSIind});

SSPind = [7,20,44];
SSP = cat(1,data{1}{SSPind});
SSPid = cat(1,animal{SSPind});

SABind = [8:10,21:23,45:49,59:62,65:67,69:70];
SAB = cat(1,data{1}{SABind});
SABid = cat(1,animal{SABind});

SAIind = [11:12,24:25,62:63];
SAI = cat(1,data{1}{SAIind});
SAIid = cat(1,animal{SAIind});

SAPind = [13,26,64];
SAP = cat(1,data{1}{SAPind});
SAPid = cat(1,animal{SAPind});

NSBind = [27:29,71:73];
NSB = cat(1,data{1}{NSBind});
NSBid = cat(1,animal{NSBind});

NSIind = [30:31,74:75];
NSI = cat(1,data{1}{NSIind});
NSIid = cat(1,animal{NSIind});

NSPind = [32,76];
NSP = cat(1,data{1}{NSPind});
NSPid = cat(1,animal{NSPind});

NABind = [33:35,50:55,77:79,83:85];
NAB = cat(1,data{1}{NABind});
NABid = cat(1,animal{NABind});

NAIind = [36:37,56:57,80:81];
NAI = cat(1,data{1}{NAIind});
NAIid = cat(1,animal{NAIind});

NAPind = [38,48,82];
NAP = cat(1,data{1}{NAPind});
NAPid = cat(1,animal{NAPind});
%% NVHL base vs. sham base
% Combine all nvhl baseline data
nvhl = [NSB;NAB];
% Combine all nvhl baseline animal IDs
nvhlID = [NSBid;NABid];
% Get list of unique nvhl IDs
uNvhl = unique(nvhlID);
% Repate above for sham baselines
sham = [SSB;SAB];
shamID = [SSBid;SABid];
uSham = unique(shamID);
% Count number of samples per animal (looks like 600 samples per will work)
shamN = cellfun(@(x) numel(logicFind(x,shamID,'==')),uSham);
nvhlN = cellfun(@(x) numel(logicFind(x,nvhlID,'==')),uNvhl);
% Determine fewest number of animals in the two data sets
n = min(numel(uNvhl),numel(uSham));
% Randomally select n animals from each data set
theseNvhlIDs = uNvhl(randperm(numel(uNvhl),n));
theseShamIDs = uSham(randperm(numel(uSham),n));
% Subset data to the same number of samples (600) per animal
% for ii = 1:n
%     this = sham(logicFind(theseShamIDs{ii},shamID,'=='),:);
%     thisSham{ii}  = this(randperm(size(this,1),600),:);
%     
%     this = nvhl(logicFind(theseNvhlIDs{ii},nvhlID,'=='),:);
%     thisNvhl{ii} = this(randperm(size(this,1),600),:);
% end
% Create matrix of all possible pairs of one sham and one nvhl
c = 1;
for ii = 1:n
    for jj = 1:n
        cmb(c,:) = [ii,jj];
        c = c+1;
    end
end
% save('D:/emily/rip/nvhlBase_shamBaseRel.mat','thisSham','thisNvhl','cmb')
%%
samps = 600;
inds = [55:78,91:96,103:108,127:132,139:162,169:174,181:204];
% Remove outlier 8th animal of sham
% uSham = uSham([1:7,9,10]);
for m = 1:100
    theseNvhlIDs = uNvhl(randperm(numel(uNvhl),n));
    theseShamIDs = uSham(randperm(numel(uSham),n));
    for ii = 1:n
        this = sham(logicFind(theseShamIDs{ii},shamID,'=='),:);
        thisSham{ii}  = this(randperm(size(this,1),samps),:);
        
        this = nvhl(logicFind(theseNvhlIDs{ii},nvhlID,'=='),:);
        thisNvhl{ii} = this(randperm(size(this,1),samps),:);
    end
    [trainX,trainY,testX,testY] = trainTest(cat(1,thisNvhl{:},...
        thisSham{:}),[ones(samps*n,1);zeros(samps*n,1)],0.2);
    if m <50
        thisCmb = cmb(m,:);
        nvhlTrain = ~ismember(1:7,thisCmb(1));
        shamTrain = ~ismember(1:7,thisCmb(2));
        testX = cat(1,thisNvhl{thisCmb(1)},thisSham{thisCmb(2)});
        testY = [ones(samps,1);zeros(samps,1)];
        trainX = cat(1,thisNvhl{nvhlTrain},thisSham{shamTrain});
        trainY = [ones(samps*(n-1),1);zeros(samps*(n-1),1)];
        mdl = fitglm(zscore(trainX(:,inds)),trainY,'distribution','binomial',...
            'binomialSize',numel(trainY));
        pred = predict(mdl,zscore(testX(:,inds)));
        [~,~,~,aLOO(m)] = perfcurve(testY,pred,1);
        % Single feature
        for f = 1:size(trainX,2)
            mdl = fitglm(zscore(trainX(:,inds(f))),trainY,'distribution',...
                'binomial','binomialSize',numel(trainY));
            pred = predict(mdl,zscore(testX(:,inds(f))));
            [~,~,~,aLOOSingle(f,m)] = percurve(testY,pred,1);
        end
    end
    mdl = fitglm(zscore(trainX(:,inds)),trainY,'distribution','binomial',...
        'binomialSize',numel(trainY));
    pred = predict(mdl,zscore(testX(:,inds)));
    [~,~,~,aReal(m)] = perfcurve(testY,pred,1);
    % Single feature
    for f = 1:size(trainX,2)
        mdl = fitglm(zscore(trainX(:,inds(f))),trainY,'distribution',...
            'binomial','binomialSize',numel(trainY));
        pred = predict(mdl,zscore(testX(:,inds(f))));
        [~,~,~,aRealSingle(f,m)] = percurve(testY,pred,1);
    end
    % Animal detector
    all = cat(2,thisSham,thisNvhl);
    a = randperm(7,7);
    b = randperm(7,7)+7;
    if m <=50
        one = cat(1,all{a(1:3)},all{b(1:4)});
        zero = cat(1,all{a(4:7)},all{b(5:7)});
    else
        one = cat(1,all{a(1:4)},all{b(1:3)});
        zero = cat(1,all{a(5:7)},all{b(4:7)});
    end
    [trainX,trainY,testX,testY] = trainTest([one;zero],[ones(samps*7,1);...
        zeros(samps*7,1)],0.2);
    mdl = fitglm(zscore(trainX(:,inds)),trainY,'distribution','binomial',...
        'binomialSize',numel(trainY));
    pred = predict(mdl,zscore(testX(:,inds)));
    [~,~,~,aAD(m)] = perfcurve(testY,pred,1);
end
%%
figure
hold on
h(1) = histogram(aReal);
h(2) = histogram(aAD);
h(3) = histogram(aLOO);
xax = get(gca,'xlim');
len = xax(2)-xax(1);
% Split into 40 bins
bin = len/40;
% Apply
set(h(1),'BinWidth',bin);
set(h(2),'BinWidth',bin);
set(h(3),'BinWidth',bin);
legend({['80/20: ',num2str(round(mean(aReal),2)),'\pm',...
    num2str(round(conf(aReal,0.95),2))],...
    ['AD: ',num2str(round(mean(aAD),2)),'\pm',...
    num2str(round(conf(aAD,0.95),2))]...
    ['LOO: ',num2str(round(mean(aLOO),2)),'\pm',...
    num2str(round(conf(aLOO,0.95),2))]})
title('logistic')
%% Load nvhl base vs. sham base model
% cd 'D:/emily/rip/nvhlBase_shamBaseRel/'
% for ii = [1:29,33:49]
%     load(['nvhlBase_shamBase',num2str(ii),'.mat'],'acc','accR')
%     allA(ii) = acc{1}.acc;
%     allAR(ii) = accR{1}.acc;
% end
% doubleHist(allA,allAR)
%% combine data
% for ii = 1:100
%    load(['D:\emily\rip\new\nvhlBase_shamBase',num2str(ii),'.mat'])
%    save(['D:\emily\rip\nvhlBase_shamBase\nvhlBase_shamBase',num2str(ii),'.mat'],...
%        'lam','beta','fits','acc','hist','lamAD','betaAD','fitsAD',...
%        'accAD','histAD','-append')
% end
%%
cd 'D:/emily/rip/nvhlBase_shamBase/'
for ii = 1:100
    disp(ii)
    if ii >= 50 
        load(['nvhlBase_shamBase',num2str(ii),'.mat'],'acc','accR',...
            'accAD','hist','histAD')
    else 
        load(['nvhlBase_shamBase',num2str(ii),'.mat'],'acc','accR',...
            'accAD','accLOO','hist','histAD')
    aLOOLasso(ii) = accLOO{1}.acc;
    end
    aLasso(ii) = acc{1}.acc;
    aRLasso(ii) = accR{1}.acc;
    aADLasso(ii) = accAD{1}.acc;
    % Single feature
    for f = 1:size(hist.trainX,2)
        mdl = fitglm(zscore(hist.trainX(:,f)),hist.trainY,...
            'distribution','binomial','binomialSize',numel(hist.trainY));
        pred = predict(mdl,zscore(hist.cfg.naive.testX(:,f)));
        [~,~,~,aLassoSingle(ii,f)] = perfcurve(hist.cfg.naive.testY,...
            pred,1);
        
        mdl = fitglm(zscore(histAD.trainX(:,f)),histAD.trainY,...
            'distribution','binomial','binomialSize',numel(histAD.trainY));
        pred = predict(mdl,zscore(histAD.cfg.naive.testX(:,f)));
        [~,~,~,aADLassoSingle(ii,f)] = perfcurve(histAD.cfg.naive.testY,...
            pred,1);
    end
end
% save('D:/emily/rip/nvhlBase_shamBase_aucs.mat','aLOOLasso','aRLasso',...
%     'aADLasso','aLassoSingle','aADLassoSingle')
%%
figure
hold on
h(1) = histogram(aLasso);
h(2) = histogram(aADLasso);
h(3) = histogram(aLOOLasso);
xax = get(gca,'xlim');
len = xax(2)-xax(1);
% Split into 40 bins
bin = len/40;
% Apply
set(h(1),'BinWidth',bin);
set(h(2),'BinWidth',bin);
set(h(3),'BinWidth',bin);
legend({['80/20: ',num2str(round(mean(aLasso),2)),'\pm',...
    num2str(round(conf(aLasso,0.95),3))],...
    ['AD: ',num2str(round(mean(aADLasso),2)),'\pm',...
    num2str(round(conf(aADLasso,0.95),2))]...
    ['LOO: ',num2str(round(mean(aLOOLasso),2)),'\pm',...
    num2str(round(conf(aLOOLasso,0.95),2))]})
title('lasso')
%% single features
nameVect = names({'PLL','M2L','CCL','NAcL','CCR','NAcR','PLR','M2R'},...
    {'d','t','a','b','lg','hg'});
inds = [55:78,91:96,103:108,127:132,139:162,169:174,181:204];
nameVect = nameVect(inds);
mAD = mean(aADLassoSingle,1);
mA = mean(aLassoSingle,1);
[smA,sortInd] = sort(mA','descend');
figure
plot(mA,mAD,'.')
%% PCA
[coeff,score,latent,~,ex] = pca([nvhl;sham]);
xCentered = score*coeff';
% figure
% hold on
% scatter3(score(1:6060,1),score(1:6060,2),score(1:6060,3),'.')
% scatter3(score(6060:end,1),score(6060:end,2),score(6060:end,3),'.')
% xlabel('1'); ylabel('2'); zlabel('3')
% legend({'nvhl','sham'},'location','ne')
figure
scatterhist(score(:,1),score(:,2),'group',[repmat('nvhl',6060,1);...
    repmat('sham',10419,1)],'marker','.','Kernel','on','bandwidth',...
    [5,5;5,5],'color','gy')
hold on
for ii = 1:numel(uSham)
these = logicFind(uSham{ii},shamID,'==');
scatter3(mean(score(these,1)),mean(score(these,2)),mean(score(these,3)),2000,'.r')%c(ii,:),'.')
end
for ii = 1:numel(uNvhl)
these = logicFind(uNvhl{ii},nvhlID,'==');
scatter3(mean(score(these,1)),mean(score(these,2)),mean(score(these,3)),2000,'.b')%c(ii+10,:),'.')
end
title('all data')
%%
figure
scatterhist([nvhl(:,inds(17));sham(:,inds(17))],[nvhl(:,inds(73));sham(:,inds(73))],'group',...
    [repmat('nvhl',6060,1);repmat('sham',10419,1)],'marker','.')
%%
data = score;%[nvhl(:,inds(sortInd));sham(:,inds(sortInd))];
[countN,centerN] = hist3(score(1:6060,1:2),[100 100]);
[countS,centerS] = hist3(score(6061:10419,1:2),[100 100]);
figure
hold on
contour(centerN{1},centerN{2},countN,[3 3],'b')
contour(centerS{1},centerS{2},countS,[3 3],'r')
%%
figure
hold on
scatter(data(:,1),data(:,2),'.k')
for ii = 1:numel(uSham)
these = logicFind(uSham{ii},shamID,'==');
scatter(mean(data(these,1)),mean(data(these,2)),2000,'.r')%c(ii,:),'.')
end
for ii = 1:numel(uNvhl)
these = logicFind(uNvhl{ii},nvhlID,'==');
scatter(mean(data(these,1)),mean(data(these,2)),2000,'.b')%c(ii+10,:),'.')
end
%%
c = distinguishable_colors(17);
figure
hold on
for ii = 1:numel(uSham)
    these = logicFind(uSham{ii},shamID,'==');
    scatter3(score(these,1),score(these,2),score(these,3),[],c(ii,:),'.')
%     scatter3(mean(score(these,1)),mean(score(these,2)),mean(score(these,3)),20,c(ii,:),'.')
end
for ii = 1:numel(uSham)
these = logicFind(uSham{ii},shamID,'==');
scatter3(mean(score(these,1)),mean(score(these,2)),mean(score(these,3)),2000,'.r')%c(ii,:),'.')
end
for ii = 1:numel(uNvhl)
these = logicFind(uNvhl{ii},nvhlID,'==');
scatter3(mean(score(these,1)),mean(score(these,2)),mean(score(these,3)),2000,'.b')%c(ii+10,:),'.')
end
%%
figure
hold on
for ii = 1:numel(uSham)
these = logicFind(uSham{ii},shamID,'==');
scatter(mean(score(these,1)),mean(score(these,2)),2000,'.r')%c(ii,:),'.')
end
for ii = 1:numel(uNvhl)
these = logicFind(uNvhl{ii},nvhlID,'==');
scatter(mean(score(these,1)),mean(score(these,2)),2000,'.b')%c(ii+10,:),'.')
end
%% highlight 7 random shams
figure
hold on
chosen = randperm(10,7);
for ii = chosen
    these = logicFind(uSham{ii},shamID,'==');
    scatter3(score(these,1),score(these,2),score(these,3),[],'.k')
end
nums = 1:10;
% other = nums(~ismember(nums,chosen));
% for ii = other
%         these = logicFind(uSham{ii},shamID,'==');
%     scatter3(score(these,1),score(these,2),score(these,3),[],[0.5 0.5 0.5],'.')
% end

scatter3(score(1:6060,1),score(1:6060,2),score(1:6060,3),'.r')
