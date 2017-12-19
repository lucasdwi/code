
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat')
base = all;
inds = 1:60;
inds = inds(~ismember(inds,pInds));
cmbs = nchoosek(1:numel(inds),3);
for ii = [3:7,10,13,15,16]
    load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\dep24_500Each6000All50-50.mat')
    dep24X = all.testX{ii}(:,inds);
    dep24Y = all.testY{ii};
    load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\dep48_500Each6000All50-50.mat')
    dep48X = all.testX{ii}(:,inds);
    dep48Y = all.testY{ii};
    load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\chow500Each6000All50-50.mat')
    chowX = all.testX{ii}(:,inds);
    chowY = all.testY{ii};
    
    % Build models at baseline using those numbers of features
    trainX = base.trainX{ii}(:,inds);
    trainY = base.trainY{ii};
    testX = base.testX{ii}(:,inds);
    testY = base.testY{ii};
    for k = 1:length(cmbs)
        mdl = fitglm(trainX(:,cmbs(k,:)),trainY,'distribution','binomial');
        prob = predict(mdl,testX(:,cmbs(k,:)));
        [~,~,~,A(ii,k)] = perfcurve(testY,prob,1);
        % Test those models across conditions
        % Dep 24
        prob = predict(mdl,dep24X(:,cmbs(k,:)));
        [~,~,~,dep24A(ii,k)] = perfcurve(dep24Y,prob,1);
        % Dep 48
        prob = predict(mdl,dep48X(:,cmbs(k,:)));
        [~,~,~,dep48A(ii,k)] = perfcurve(dep48Y,prob,1);
        % Chow
        prob = predict(mdl,chowX(:,cmbs(k,:)));
        [~,~,~,chowA(ii,k)] = perfcurve(chowY,prob,1);
        if sum(k == round((0.05:0.05:1).*30856)) >= 1
            disp([num2str(ii),'-',num2str(k/30856*100),'%'])
        end
    end
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\trip',num2str(ii),'.mat'],'A','dep24A','dep48A','chowA')
end