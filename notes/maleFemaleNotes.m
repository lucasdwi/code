%% Collate all trls
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\maleFemale\processed\'],{'estrus';'di';'pro';'male'},...
    {'pow','coh'},'trl','rel');
% Combine each into indpendent matrices
est = cat(1,data{1}{:});
di = cat(1,data{2}{:});
pro = cat(1,data{3}{:});
male = cat(1,data{4}{:});
% Re-combine into one cell array for easy iterations
catData{1} = est;
catData{2} = di;
catData{3} = pro;
catData{4} = male;
%% Build models for each pair of cells
cmbs = nchoosek(1:4,2);
samp = cellfun(@(x) size(x,1),catData);
for ii = 1:size(cmbs,1)
   % Determine minimum samples
   minSamp = min([samp(cmbs(ii,1)),samp(cmbs(ii,2))]);
   % Determine 80/20 split based on minSamp
   trainN = round(minSamp*.8); testN = round(minSamp*.2);
   for jj = 1:100
       if rem(jj,10) == 0
           disp([num2str(ii),'.',num2str(jj)])
       end
       % Get indices for both sets
       trainInd = [randperm(samp(cmbs(ii,1)),trainN);...
           randperm(samp(cmbs(ii,2)),trainN)];
       testInd = [randperm(samp(cmbs(ii,1)),testN);...
           randperm(samp(cmbs(ii,2)),testN)];
       % Prep train and test x/y
       trainX = [catData{cmbs(ii,1)}(trainInd(1,:),:);...
           catData{cmbs(ii,2)}(trainInd(2,:),:)];
       trainY = [zeros(trainN,1);ones(trainN,1)];
       testX = [catData{cmbs(ii,1)}(testInd(1,:),:);...
           catData{cmbs(ii,2)}(testInd(2,:),:)];
       testY = [zeros(testN,1);ones(testN,1)];
       % Build model
       mdl{ii,jj} = fitglm(trainX,trainY,'distribution','binomial');
       prob = predict(mdl{ii,jj},testX);
       [x{ii}(jj,:),y{ii}(jj,:),~,a(ii,jj)] = perfcurve(testY,prob,1);
       % Also test on permuted data
       [xRand{ii}(jj,:),yRand{ii}(jj,:),~,aRand(ii,jj)] = perfcurve(...
           testY(randperm(size(testY,1))),prob,1);
       for k = 1:60
           % Build model
           logMdl{ii,jj,k} = fitglm(trainX(:,k),trainY,'distribution',...
               'binomial');
           prob = predict(logMdl{ii,jj,k},testX(:,k));
           [logX{ii,k}(jj,:),logY{ii,k}(jj,:),~,logA(ii,jj,k)] = ...
               perfcurve(testY,prob,1);
           % Also test on permuted data
           [logXRand{ii,k}(jj,:),logYRand{ii,k}(jj,:),~,logARand(ii,jj,k)] = ...
               perfcurve(testY(randperm(size(testY,1))),prob,1);
       end
   end
end
%%
names = {'est','di','pro','male'};
figure
for ii = 1:6
    subplot(2,3,ii)
    hold on
    plot(mean(x{ii},1),mean(y{ii},1))
    plot(mean(xRand{ii},1),mean(yRand{ii},1),'--k')
    title([names{cmbs(ii,1)},' vs. ',names{cmbs(ii,2)}])
    legend({['Actual: ',num2str(round(mean(a(ii,:)),2)),'\pm',...
        num2str(round(conf(a(ii,:),0.95),4))],...
        ['Permuted: ',num2str(round(mean(aRand(ii,:)),2)),'\pm',...
        num2str(round(conf(aRand(ii,:),0.95),4))]},'location','se')
end
%% Estrus vs. Diestrus
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\maleFemale\processed\'],{'estrus';'di'},{'pow','coh'},'avg','');
estrus = cat(1,data{1,1}{:});
diestrus = cat(1,data{1,2}{:});
data = [estrus;diestrus];
y = [ones(size(estrus,1),1);zeros(size(diestrus,1),1)];
%% Estrus vs. Male
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'estrus';'male'},{'pow','coh'},'avg','');
estrus = cat(1,data{1,1}{:});
male = cat(1,data{1,2}{:});
data = [estrus;male];
y = [ones(size(estrus,1),1);zeros(size(male,1),1)];
% for ii = 1:size(data,2)
%     mdl = fitglm(data(:,ii),y,'distribution','binomial');
%     r(ii) = mdl.Rsquared.Ordinary;
% end
%% Female vs. Male
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'di';'male'},{'pow','coh'},'avg','');
female = cat(1,data{1,1}{:});
male = cat(1,data{1,2}{:});
data = [female;male];
y = [ones(size(female,1),1);zeros(size(male,1),1)];
%% Estrus vs. Pro
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'estrus';'pro'},{'pow','coh'},'avg','');
estrus = cat(1,data{1,1}{:});
pro = cat(1,data{1,2}{:});
data = [estrus;pro];
y = [ones(size(estrus,1),1);zeros(size(pro,1),1)];
%% Diestrus vs. Pro
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'di';'pro'},{'pow','coh'},'avg','');
diestrus = cat(1,data{1,1}{:});
pro = cat(1,data{1,2}{:});
data = [diestrus;pro];
y = [ones(size(diestrus,1),1);zeros(size(pro,1),1)];
%% Pro vs. Male
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'pro';'male'},{'pow','coh'},'avg','');
pro = cat(1,data{1,1}{:});
male = cat(1,data{1,2}{:});
data = [pro;male];
y = [ones(size(pro,1),1);zeros(size(male,1),1)];
%% Load data by trials - for getting first 10 minutes of pre-EtOH data
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkNot\',{'.mat'},{'pow','coh'},'trl','rel');
% Load alcohol amounts
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\fileAlc.mat')
% Go through trials and only keep those occuring within first 10 minutes
% Set number of minutes to collect up to
preTime = 10;
% Convert to samples
preSamp = preTime*60*400;
% Preallocate
preData = zeros(size(data{1},1),size(data{1}{1,2},2));
for ii = 1:size(data{1},1)
    % Find index of last trial before preSamp
    ind = logicFind(preSamp,samp{1}{ii,2}(:,2),'<=','last');
    % Average those trials together and concatenate to preData
    preData(ii,:) = mean(data{1}{ii,2}(1:ind,:),1);
    % Split name of current file
    splits = strsplit(files{1}{ii},'_d');
    % Go through files and find corresponding drinking amount
    ind = find(not(cellfun(@isempty,strfind(alc(:,1),splits{1}))));
    % Grab alcohol amount
    yAlc(ii,1) = alc{ind,2};
end
% save('C:\Users\Pythia\Documents\GreenLab\data\paper3\preDataAlc.mat','preData','yAlc')
%% Build model with only males and only female (diestrus); test across
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\preDataAlc.mat')
group = cellstr([repmat('m',18,1);'d';'p';'e';'d';'d';'m';'m';'d';'p';...
    'd';'e';'p';'d';'d';'d';'d';'e';'d';'d';'d';'d';'p';'d';'p';'e';'d';...
    'p']);
mInd = logicFind('m',group,'==');
dInd = logicFind('d',group,'==');
maleData = preData(mInd,:);
maleY = yAlc(mInd);
diData = preData(dInd,:);
diY = yAlc(dInd);
save('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\maleDiAlc.mat',...
    'maleData','maleY','diData','diY')
%% Predicting drinking amounts - all files (male and all female)
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'preDrinkAlc.mat'])
real = allLambda{1}.allErr;
perm = cellfun(@(x) x.allErr, allLambdaPerm,'UniformOutput',0);
perm = cat(1,perm{:});
doubleHist(real,perm,'Main','Predicting alcohol consumption','xlab',...
    'MAE (g/kg)');
%% Predicting drinking amounts - males
load(['C:\Users\Pythia\Documents\GreenLab\data\maleFemale\'...
    'preDrinkAlcMale.mat'])
real = allLambda{1}.allErr;
perm = cellfun(@(x) x.allErr,allLambdaPerm,'UniformOutput',0);
perm = cat(1,perm{:});
doubleHist(real,perm,'Main','Predicting alcohol consumption: Males',...
    'xlab','MAE (g/kg)');
%% Predicting drinking amounts - diestrus female
load(['C:\Users\Pythia\Documents\GreenLab\data\maleFemale\'...
    'preDrinkAlcDi.mat'])
real = allLambda{1}.allErr;
perm = cellfun(@(x) x.allErr,allLambdaPerm,'UniformOutput',0);
perm = cat(1,perm{:});
doubleHist(real,perm,'Main','Predicting alcohol consumption: Diestrus',...
    'xlab','MAE (g/kg)');

%%
load('C:\Users\Lucas\Desktop\maleFemale\phase1.mat')
realAcc = (1-real{1}.err).*100;
rM = mean(realAcc);
rS = std(realAcc);
permAcc = (1-perm{1}.err).*100;
pM = mean(permAcc);
pS = std(permAcc);
[d] = distES(realAcc,permAcc);
%%
figure
histogram(realAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','k','EdgeColor','w')
hold on
histogram(permAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','w')
yticks = get(gca,'Ytick').*100;
set(gca,'YTickLabel',yticks);
ylabel('Percent of Models (%)')
xlabel('Accuracy')
title('Diestrus vs. Estrus')
legend({['Real: ',num2str(round(rM)),'\pm',num2str(round(rS)),'%'],['Permuted: ',num2str(round(pM)),'\pm',num2str(round(pS)),'%']},'Location','northwest')
text(22,.086,['d = ',num2str(round(d,2))])
box off
%%
load('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\maleFemale1.mat')
realAcc = (1-real{1}.err).*100;
rM = mean(realAcc);
rS = std(realAcc);
permAcc = (1-perm{1}.err).*100;
pM = mean(permAcc);
pS = std(permAcc);
[d] = distES(realAcc,permAcc);
%%
figure
histogram(realAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','k','EdgeColor','w')
hold on
histogram(permAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','w')
title('Male vs. Female (Estrus)')
ylabel('Percent of Models (%)')
yticks = get(gca,'Ytick').*100;
set(gca,'YTickLabel',yticks);
xlabel('Accuracy (%)')
legend({['Real: ',num2str(round(rM)),'\pm',num2str(round(rS)),'%'],['Permuted: ',num2str(round(pM)),'\pm',num2str(round(pS)),'%']},'Location','northwest')
text(12,.076,['d = ',num2str(round(d,2))])
box off
%%
load('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\maleDi1.mat')
realAcc = (1-real{1}.err).*100;
rM = mean(realAcc);
rS = std(realAcc);
permAcc = (1-perm{1}.err).*100;
pM = mean(permAcc);
pS = std(permAcc);
[d] = distES(realAcc,permAcc);
%%
figure
histogram(realAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','k','EdgeColor','w')
hold on
histogram(permAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','w')
title('Male vs. Female (Diestrus)')
ylabel('Percent of Models (%)')
yticks = get(gca,'Ytick').*100;
set(gca,'YTickLabel',yticks);
xlabel('Accuracy (%)')
legend({['Real: ',num2str(round(rM)),'\pm',num2str(round(rS)),'%'],['Permuted: ',num2str(round(pM)),'\pm',num2str(round(pS)),'%']},'Location','northwest')
text(22,.086,['d = ',num2str(round(d,2))])
box off
%%
load('maleDi.mat')
cmbs = nchoosek(1:26,5);
inds = 1:26;
for ii = 1:100
    chk = 1;
    while chk == 1
        samp = cmbs(randi(65780,1,1),:);
        chk = size(unique(y(samp)),1);
    end
    test(:,ii) = samp;
end

for ii = 1:size(data,2)
    for jj = 1:100
        train = inds(~ismember(inds,test(:,jj)));
        mdl = fitglm(data(train,ii),y(train),'distribution','binomial');
        prob = predict(mdl,data(test(:,jj),ii));
        dir(ii,jj) = exp(table2array(mdl.Coefficients(2,1)));
        [~,~,~,a(ii,jj)] = perfcurve(y(test(:,jj)),prob,1);
    end
end
%

mA = mean(a,2);
sA = std(a,[],2);
[smA,inds] = sort(mA,'descend');
tierStart = tier(a(:,inds));
ssA = sA(inds);
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
nameSort = nameVect(inds)';
%%
scatterErr(1:60,smA,ssA,1)
hold on
% Plot 50 line
plot([0 60],[0.5 0.5],'--k','LineWidth',2)

xlabel('Feature')
ylabel('AUC')
title('Average AUC from Univariate Logistic')
%%
load('maleFemale.mat')
cmbs = nchoosek(1:20,5);
inds = 1:20;
for ii = 1:100
    chk = 1;
    while chk == 1
        samp = cmbs(randi(size(cmbs,1),1,1),:);
        chk = size(unique(y(samp)),1);
    end
    test(:,ii) = samp;
end

for ii = 1:size(data,2)
    for jj = 1:100
        train = inds(~ismember(inds,test(:,jj)));
        mdl = fitglm(data(train,ii),y(train),'distribution','binomial');
        prob = predict(mdl,data(test(:,jj),ii));
        dir(ii,jj) = exp(table2array(mdl.Coefficients(2,1)));
        [~,~,~,a(ii,jj)] = perfcurve(y(test(:,jj)),prob,1);
    end
end
%%
mA = mean(a,2);
sA = std(a,[],2);
[smA,inds] = sort(mA,'descend');
ssA = sA(inds);
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
nameSort = nameVect(inds)';
%%
load('phase.mat')
cmbs = nchoosek(1:26,5);
inds = 1:26;
for ii = 1:100
    chk = 1;
    while chk == 1
        samp = cmbs(randi(65780,1,1),:);
        chk = size(unique(y(samp)),1);
    end
    test(:,ii) = samp;
end

for ii = 1:size(data,2)
    for jj = 1:100
        train = inds(~ismember(inds,test(:,jj)));
        mdl = fitglm(data(train,ii),y(train),'distribution','binomial');
        prob = predict(mdl,data(test(:,jj),ii));
        dir(ii,jj) = exp(table2array(mdl.Coefficients(2,1)));
        [~,~,~,a(ii,jj)] = perfcurve(y(test(:,jj)),prob,1);
    end
end
%%
mA = mean(a,2);
sA = std(a,[],2);
[smA,inds] = sort(mA,'descend');
ssA = sA(inds);
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
nameSort = nameVect(inds)';
scatterErr(1:60,smA,ssA,1)
hold on
% Plot 50 line
plot([0 60],[0.5 0.5],'--k','LineWidth',2)

xlabel('Feature')
ylabel('AUC')
title('Average Univariate AUC: Estrus vs. Diestrus')