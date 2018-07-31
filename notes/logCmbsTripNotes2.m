load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\concatData.mat')
cmbs = nchoosek(1:60,3);
for ii = 3:20
    disp(num2str(ii))
    tic
    trainX = all.trainX{ii};
    trainY = all.trainY{ii};
    testX = all.testX{ii};
    testY = all.testY{ii};
    for k = 1:length(cmbs)
        mdl = fitglm(trainX(:,cmbs(k,:)),trainY,'distribution','binomial');
        prob = predict(mdl,testX(:,cmbs(k,:)));
        [rocX{k},rocY{k},~,a(k)] = perfcurve(testY,prob,1);
    end
    trainX = []; trainY = []; testX = []; testY = [];
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\logCmbsTrip',num2str(ii),'.mat'],'rocX','rocY','a')
    toc
end
