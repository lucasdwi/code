load('D:/emily/f31/emilyNVHLdata.mat','NAI')
cd D:/emily/f31/nvhlSham/
for ii = 1:100
    load(['emilyNVHL',num2str(ii),'.mat'])
    allAcc(ii) = acc{1}.acc;
    allAccR(ii) = accR{1}.acc;
    allA(ii,:) = a;
    allS(ii,:) = s;
    [allX(ii,:),allY(ii,:),~,~] = perfcurve(hist.cfg.naive.testY,acc{1}.pred,1);
    
    % Replace NVHL with NAI in test set
    testY = hist.cfg.naive.testY;
    testX = hist.cfg.naive.testX;
    
    testX(testY==1,:) = [];
    testY(testY==1) = [];
        
    testY = [testY;ones(6,1)];
    testX = [testX;NAI];
    
    prob = cvglmnetPredict(acc{1}.mdl{1},testX,'lambda_1se','response');
    [thisX,thisY,~,alcA(ii)] = perfcurve(testY,prob,1);
    alcX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,17));
    alcY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,17));
end
figure
hold on
plot(mean(allX,1),mean(allY,1))
plot(mean(alcX,1),mean(alcY,1))
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
legend({['NVHL v. sham: ',num2str(round(mean(allAcc),2))],['NVHL+ALC v. sham: ',num2str(round(mean(alcA),2))]})

mAllA = mean(allA,1);
mAllS = mean(allS,1);
[sortA,sortInd] = sort(mAllA,'descend');
nameVect = names({'PLL','M2L','CCL','NAcL','CCR','NAcR','PLR','M2R'},...
    {'d','t','a','b','lg','hg'});
sortFeats = nameVect(sortInd)';

cohFeat = nameVect(sortInd(sortInd>48))';
cohA = mAllA(sortInd(sortInd>48));
cohInds = sortInd(sortInd>48);

load('D:/emily/f31/emilyNVHLdataRaw.mat')
nvhl = [NSB;NAB]; sham = [SSB;SAB];

inds = [3,12,16,21,22,24,25];
figure
hold on
c=1;
for ii = inds
    plot([c c],mean(nvhl(:,cohInds(ii))),'ro')
    plot([c c],mean(sham(:,cohInds(ii))),'bo')
    plot([c c],mean(NAI(:,cohInds(ii))),'yo')
    c=c+1;
end
set(gca,'xtick',1:c-1,'xticklabel',cohFeat(inds))
xtickangle(45)
% legend({'nvhl','sham'})
%%
inds = [61:66,73:84,91:96,103:108,121:126,133:138,151:156,163:174,193:204];
for ii = 1:numel(inds)
    [~,p(ii)] = ttest2(abs(nvhl(:,inds(ii))),abs(sham(:,inds(ii))));
end
pAdj = p*numel(p);

sigP = p<=0.05;
figure
hold on
c=1;
for ii = 1:numel(sigP)
    if sigP(ii)==1
        plot([c c],mean(abs(nvhl(:,inds(ii)))),'ro')
        plot([c c],mean(abs(sham(:,inds(ii)))),'bo')
        plot([c c],mean(abs(NAI(:,inds(ii)))),'ko')
        c = c+1;
    end
end
set(gca,'xtick',1:c-1,'xticklabel',nameVect(inds(sigP)))
xtickangle(45)
%% Get mean and std of each group
nvhlM = mean(abs(nvhl(:,inds(sigP))));
nvhlS = std(abs(nvhl(:,inds(sigP))));
nvhlSE = nvhlS./size(nvhl,1);
shamM = mean(abs(sham(:,inds(sigP))));
shamS = std(abs(sham(:,inds(sigP))));
shamSE = shamS/size(sham,1);
alcM = mean(abs(NAI(:,inds(sigP))));
alcS = std(abs(NAI(:,inds(sigP))));
alcSE = alcS/size(NAI,1);
n = [nvhlM',nvhlS',nvhlSE'];
s = [shamM',shamS',shamSE'];
a = [alcM',alcS',alcSE'];