[beta{fi},dev{fi},stats{fi}] = glmfit(T1{1,fi}(:,logicFind(0.9,masterSurvBeta{1,1},'>=')+1),T1{1,fi}(:,1),'binomial','link','logit');
for iter = 1:size(T1{1,fi},1)
    vars = [1,T1{1,fi}(iter,logicFind(0.9,masterSurvBeta{1,1},'>=')+1)];
    betax = sum(vars.*beta{52}');
    prob(fi,iter) = exp(betax)/(1+exp(betax));
    if iter <= numel(tsRest{fi})
        labels{iter} = 'NotBinge';
    else
        labels{iter} = 'Binge';
    end
end
%%
c = 1;
for ti = 0:0.01:1
    bingeT = [tsBinge{1,fi},ones(size(tsBinge{1,fi},1),1),threshRound(prob(fi,end-size(tsBinge{1,fi},1)+1:end),ti,[0,1])',prob(fi,end-size(tsBinge{1,fi},1)+1:end)'];
    restT = [tsRest{1,fi},zeros(size(tsRest{1,fi},1),1),threshRound(prob(fi,1:size(tsRest{1,fi},1)),ti,[0,1])',prob(fi,1:size(tsRest{1,fi},1))'];
    allT{c} = sortrows([bingeT;restT]);
    thisT = allT{c};
    allAcc(c) = mean(thisT(:,2)==thisT(:,3));
    bingeAcc(c) = mean(thisT(thisT(:,2)==1,2)==thisT(thisT(:,2)==1,3));
    restAcc(c) = mean(thisT(thisT(:,2)==0,2)==thisT(thisT(:,2)==0,3));
    avgBingeP(c) = mean(thisT(thisT(:,2)==1,4));
    falsePos(c) = sum(thisT(thisT(:,2)==0,2)~=thisT(thisT(:,2)==0,3))/numel(thisT(thisT(:,2)==0,2));
    truePos(c) = sum(thisT(thisT(:,2)==1,2)==thisT(thisT(:,2)==1,3))/numel(thisT(thisT(:,2)==1,2));
    c = c+1;
end
%%
% Plots at 50% threshold
figure; imagesc(allT{51}(:,2:4))
% Calculates AUC and Plots ROC
figure; plot(falsePos,truePos)
xlabel('False Pos'); ylabel('True Pos');