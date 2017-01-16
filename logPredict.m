nameVect = names;
%%
for fi = 2%:length(T1)
    [masterDev{fi},masterLam{fi},masterBeta{fi},masterMinLam{fi},masterSurvBeta{fi},masterMeanBeta{fi},masterOR{fi},masterStdBeta{fi},TOrig{fi},betaNames{fi},hist{fi}] = params(T1{1,fi},'y',1,'n','binom',[1],0.9,10,'n',(1),nameVect);
    save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\enetAutosave.mat')
end
    [beta{fi},dev{fi},stats{fi}] = glmfit(T1{1,fi}(:,logicFind(0.99,masterSurvBeta{fi}{1,1},'>=')+1),T1{1,fi}(:,1),'binomial','link','logit');
    %md1{fi} = fitglm(zscore(T1{1,2}(:,logicFind(0.99,masterSurvBeta{1,1},'>=')+1)),T1{1,2}(:,1),'distribution','binomial','link','logit')
    for iter = 1:size(T1{1,fi},1)
        vars = [1,T1{1,fi}(iter,logicFind(0.99,masterSurvBeta{1,fi}{1,1},'>=')+1)];
        betax = sum(vars.*beta{fi}');
        prob(fi,iter) = exp(betax)/(1+exp(betax));
    end
    bingeT = [tsBinge{1,fi},ones(size(tsBinge{1,fi},1),1),round(prob(end-size(tsBinge{1,fi},1)+1:end))',prob(end-size(tsBinge{1,fi},1)+1:end)'];
    restT = [tsRest{1,fi},zeros(size(tsRest{1,fi},1),1),round(prob(1:size(tsRest{1,fi},1)))',prob(1:size(tsRest{1,fi},1))'];
    allT = sortrows([bingeT;restT]);
end