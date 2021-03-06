%% New aaberg analysis with animals from Diana
[data,~,files] = collateData('E:\aaberg\processedNew\three\',...
    {'Control','in';'THC','in'},{'pow','coh'},'avg','rel');
% Combine into datasets with either 3 or two recordings per animal
% Three
x3 = cat(1,data{1}{:},data{2}{:});
y3 = [zeros(21,1);ones(30,1)];
% Two
inds0 = [1,2,4,5,7,8,10,11,13,14,16,17,19,20];
inds1 = [1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29];
x2 = cat(1,data{1}{inds0},data{2}{inds1});
y2 = [zeros(14,1);ones(20,1)];
save('E:\aaberg\newModelData.mat','x3','y3','x2','y2')
%%
[data,~,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\modelData\',{'AC_','in';'AT_','in'},{'pow','coh'},'avg','rel');
%% Open lasso results
load('E:\aaberg\newModel_5fold.mat')
doubleHist((1-real2.err).*100,(1-perm2.err).*100,'xlab','Accuracy (%)','main','2 sample')
doubleHist((1-real3.err).*100,(1-perm3.err).*100,'xlab','Accuracy','main','3 sample')
%%
for ii = 1:70
    load(['E:\aaberg\doubleDip\newModel_5fold_doub',num2str(ii),'.mat'])
    aRD(ii) = real2.acc{1}.acc;
    aPD(ii) = perm2.acc{1}.acc;
    cvD(ii,:) = real2.err;
    cvPD(ii,:) = perm2.err;
    
    load(['E:\aaberg\tripleDip\newModel_5fold_trip',num2str(ii),'.mat'])
    aRT(ii) = real3.acc{1}.acc;
    aPT(ii) = perm3.acc{1}.acc;
    cvT(ii,:) = real3.err;
    cvPT(ii,:) = perm3.err;
end
%%
load('E:\aaberg\newModelData.mat')
for ii = 1:70
    load(['E:\aaberg\tripleDipRandom\newModel_5fold_trip',num2str(ii),...
        '.mat'])
    aRT(ii) = real3.acc{1}.acc;
    aPT(ii) = perm3.acc{1}.acc;
    cvT(ii,:) = real3.err;
    cvPT(ii,:) = perm3.err;
    [x(ii,:),y(ii,:),~,a(ii)] = perfcurve(y3(real3.hist.testInd),...
        real3.acc{1}.pred,1,'TVals',0:0.05:1,'UseNearest',0);
    [xP(ii,:),yP(ii,:),~,aP(ii)] = perfcurve(y3(perm3.hist.testInd),...
        perm3.acc{1}.pred,1,'TVals',0:0.05:1,'UseNearest',0);
end
figure
plot(mean(x,1),mean(y,1),'-k')
hold on
plot(mean(xP,1),mean(yP,1),'--k')
box off
title('THC vs. Control')
xlabel('FPR'); ylabel('TPR')
legend({['Real: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),2))],['Permuted: ',...
    num2str(round(mean(aP),2)),'\pm',num2str(round(conf(aP,0.95),2))]})
%% Run single feature - leave two animals out (triple)
load('E:\aaberg\newModelData.mat')
zero = reshape(1:21,3,7)';
one = reshape(22:51,3,10)';
x = [];
for ii = 1:7
    for jj = 1:10
        x = [x;ii,jj];
    end
end
for ii = 1:size(x,1)
    disp(num2str(ii))
    testInds = [zero(x(ii,1),:),one(x(ii,2),:)];
    allInds = 1:51;
    trainInds = ~ismember(allInds,testInds);
    for jj = 1:216
        mdl = fitglm(x3(trainInds,jj),y3(trainInds),'distribution',...
            'binomial','binomialsize',sum(trainInds));
        slope(ii,jj) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,x3(testInds,jj));
        [~,~,~,a(ii,jj)] = perfcurve(y3(testInds),prob,1);
    end
end
mA = mean(a,1);
[mASort,sortInd] = sort(mA,'descend');
mASort = mASort';
slopeSort = mean(slope(:,sortInd),1)';
nameVect = names({'rOFC','rIL','rNAc','rPL','lNAc','lPL','lOFC','lIL'},...
    {'d','t','a','b','lg','hg'});
feat = nameVect(sortInd)';
%%
x = [data{1};data{2}];
y = [zeros(size(data{1},1),1);ones(size(data{2},1),1)];
allInds = 1:size(x,1);
for ii = 1:size(x,1)/3
    testInds = (3*ii-2):(3*ii);
    trainInds = allInds(~ismember(1:size(x,1),testInds));
    looTrainX{ii} = cat(1,x{trainInds,:});
    looTrainY{ii} = y(trainInds,:);
    looTestX{ii} = cat(1,x{testInds,:});
    looTestY{ii} = y(testInds,:);
end
%%
inds = [];
for ii = 1:4
    for jj = 5:9
        inds = [inds;ii,jj];
    end
end
allInds = 1:size(x,1);
for ii = 1:20
    out = [inds(ii,1)*3-2:inds(ii,1)*3,inds(ii,2)*3-2:inds(ii,2)*3];
    in = allInds(~ismember(allInds,out));
    ltoTrainX{ii} = cat(1,x{in,:});
    ltoTrainY{ii} = y(in,:);
    ltoTestX{ii} = cat(1,x{out,:});
    ltoTestY{ii} = y(out,:);
end
save('C:\Users\Pythia\Documents\GreenLab\data\aaberg\aabergModelDataLto.mat','ltoTrainX','ltoTrainY','ltoTestX','ltoTestY')
%%
for ii = 1:20
   load(['C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\lto\lto',num2str(ii),'.mat']) 
   A(ii) = a;
end
%%
load('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\allDataModel\1.mat')
doubleHist((1-real.err).*100,(1-perm.err).*100,'Main','THC vs. CON','xlab','Accuracy (%)')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\log\log1.mat')
r = cellfun(@(x) x.Rsquared.Adjusted,mdl);
nameVect = names({'rOFC','rIL','rNAc','rPL','lNAc','lPL','lOFC','lIL'},{'d','t','a','b','lg','hg'});
[sR,sInd] = sort(r);
nameSort = nameVect(sInd);
%%
load('C:\Users\Pythia\Documents\GreenLab\data\aaberg\leverModel.mat')
for ii = 3%1:4
   figure
   mdl = fitglm(x(1:3:27,sInd(end-(ii-1))),lever(1:3:27));
   scatter(x(1:3:27,sInd(end-(ii-1))),lever(1:3:27),'ok','filled')
   lsline
   title([nameVect(sInd(end-(ii-1)))])
   xAx = get(gca,'XTick');
   text(xAx(end-1),55,['R^2 = ',num2str(round(mdl.Rsquared.Ordinary,2)),newline,'p = ',num2str(round(table2array(mdl.Coefficients(2,4)),2))])
   xlabel('Normalized Coherence')
   ylabel('Average Lever Presses')
end
%% Grab all power and coherence
files{1} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\modelData\','AC_','in');
files{2} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\modelData\','AT_','in');
p = cell(1,2);
c = cell(1,2);
for ii = 1:size(files,2)
   for jj = 1:size(files{ii},2)
       load(files{ii}{jj})
       p{ii} = cat(3,p{ii},psdTrls{1,1}.Overall);
       c{ii} = cat(3,c{ii},mean(coh{1,1}.Cxy,3));
   end
end
% Average
mP = cellfun(@(x) mean(x,3),p,'UniformOutput',0);
mC = cellfun(@(x) mean(x,3),c,'UniformOutput',0);
% STD
sP = cellfun(@(x) std(x,[],3),p,'UniformOutput',0);
sC = cellfun(@(x) std(x,[],3),c,'UniformOutput',0);
% Grab frequency axes
cF = coh{1}.f;
pF = psdTrls{1}.f;
%% Plot rPL-lOFC coherence- with area under curve
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    inds = [nearest_idx3(bands{ii,2}(1),cF) nearest_idx3(bands{ii,2}(2),cF)];
    h = area(cF(inds(1):inds(2)),mC{2}(21,inds(1):inds(2)));
    h.FaceColor = [0.85 0.85 0.85];
    
    h = area(cF(inds(1):inds(2)),mC{1}(21,inds(1):inds(2)));
    h.FaceColor = [0.75 0.75 0.75];
    
    text(sum(bands{ii,2})/2,0.03,latinToGreek({bands{ii,1}}),'HorizontalAlignment','center')
end
clear h
h(1) = plot(cF,mC{1}(21,:),'-r');
h(2) = plot(cF,mC{2}(21,:),'-b');
legend(h(1:2),{'Control','THC'})

title('lNAc-lOFC Coherence')
xlabel('Frequency (Hz)')
ylabel('Magnitude Squared Coherence')
ylim([0 1])
%% Plot lIL power - with area under curve
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    inds = [nearest_idx3(bands{ii,2}(1),pF) nearest_idx3(bands{ii,2}(2),pF)];
    h = area(pF(inds(1):inds(2)),mP{2}(8,inds(1):inds(2)));
    h.FaceColor = [0.85 0.85 0.85];
    
    h = area(pF(inds(1):inds(2)),mP{1}(8,inds(1):inds(2)));
    h.FaceColor = [0.75 0.75 0.75];
    
    text(sum(bands{ii,2})/2,-3,latinToGreek({bands{ii,1}}),'HorizontalAlignment','center')
end
clear h
h(1) = plot(pF,mP{1}(8,:),'-r');
h(2) = plot(pF,mP{2}(8,:),'-b');
legend(h,{'Control','THC'},'Location','sw')
title('lIL Power')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
%% Plot rPL-lOFC coherence - with rectangle
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    rectangle('Position',[bands{ii,2}(1) 0.3 diff(bands{ii,2}) 0.47],'FaceColor',[0.85 0.85 0.85])
end
plot(cF,mC{1}(21,:),'-r')
plot(cF,mC{2}(21,:),'-b')
legend({'Control','THC'})
title('lNAc-lOFC Coherence')
xlabel('Frequency (Hz)')
ylabel('Magnitude Squared Coherence')
ylim([0 1])

%% Plot lIL power - with rectangle
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    rectangle('Position',[bands{ii,2}(1) -85 diff(bands{ii,2}) 50],'FaceColor',[0.85 0.85 0.85])
end
plot(pF,mP{1}(8,:),'-b')
plot(pF,mP{2}(8,:),'-r')
legend({'Control','THC'})
title('lIL Power')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')