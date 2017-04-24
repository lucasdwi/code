badInds = [35,48,56];
goodInds = ~ismember(1:60,badInds);
for ii = 1:size(allData,1)
    % First build model and test on self
    cfg = lassoNetCfg(0.20,'n','y','n',100,'1se');
    thisDat = allData{ii,1}(:,goodInds);
    [~,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(thisDat,resp{ii,1},'binomial','auc',1,5,1,cfg);
    % Find internal accuracy
    [selfX{ii},selfY{ii},~,selfA{ii}] = perfcurve(resp{ii,1}(hist.testInd),accArray{1,1}.pred,1);
    % Find accuracy when applied to other animals
    for oi = 1:size(allData,1)
        thisDat = allData{oi,1}(:,goodInds);
        [predY] = cvglmnetPredict(cvFitsArray{1,1}{allLambda{1,1}.bestLambdaInds},thisDat,'lambda_1se','response');
        [otherX{ii,oi},otherY{ii,oi},~,otherA(ii,oi)] = perfcurve(resp{oi,1},predY,1);
    end
end
%% Test model built in animal and tested on different conditions within animal
badInds = [35,48,56];
goodInds = ~ismember(1:60,badInds);
for ii = 1:size(allData,1)
    for j = 1:sum(~cellfun(@isempty,allData(ii,:)))
        % First build model and test on self
        cfg = lassoNetCfg(0.20,'n','y','n',100,'1se');
        thisDat = allData{ii,j}(:,goodInds);
        [~,allLambda,allBeta{ii,j},cvFitsArray,accArray,hist] = lassoNet(thisDat,resp{ii,j},'binomial','auc',1,5,1,cfg);
        models{ii,j} = cvFitsArray{1,1}{allLambda{1,1}.bestLambdaInds};
        % Find internal accuracy
        [selfX{ii,j},selfY{ii,j},~,selfA{ii,j}] = perfcurve(resp{ii,j}(hist.testInd),accArray{1,1}.pred,1);
        % Find accuracy when applied to other conditions
        for oi = 1:sum(~cellfun(@isempty,allData(ii,:)))
            if oi ~= j
                testDat = allData{ii,oi}(:,goodInds);
                [predY{ii,oi}] = cvglmnetPredict(models{ii,j},testDat,'lambda_1se','response');
                [allX{ii,oi},allY{ii,oi},~,allA(ii,oi)] = perfcurve(resp{ii,oi},predY{ii,oi},1);
            else
                allX{ii,oi} = selfX{ii,oi};
                allY{ii,oi} = selfY{ii,oi};
                allA(ii,oi) = selfA{ii,oi};
            end
        end
    end
end
%% Replace diagonal
allAUC = otherA;
for ii = 1:12
    allAUC(ii,ii) = selfA{1,ii};
end
%%
figure
pcolor([allAUC,zeros(12,1);zeros(1,13)])
xlabel('Test Set')
ylabel('Train Set')
colormap('viridis')
%%
figure
hold on
for ii = 1:12
   if ii == 1
       plot(selfX{ii},selfY{ii})
   else
        plot(otherX{1,ii},otherY{1,ii})
   end
end