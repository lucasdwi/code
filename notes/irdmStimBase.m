function runIRDMstimBase(n)
addpath(genpath('~/data'))
addpath(genpath('~/code'))
% Load data
load('~/data/stimBase.mat')
%%
for ii = 1:size(inter,2)
    for jj = 1:size(stim,1)
        thisStim{ii,jj} = cat(1,stim{jj,ii,:});
        thisBase{ii,jj} = [];
        for k = 1:size(base,3)
            if ~isempty(base{jj,ii,k}) 
                if isempty(thisBase{ii,jj})
                   thisBase{ii,jj} = [thisBase{ii,jj};base{jj,ii,k}];
                elseif base{jj,ii,k}(end) ~= thisBase{ii,jj}(end)
                    thisBase{ii,jj} = [thisBase{ii,jj};base{jj,ii,k}];
                end                    
            end
        end
    end
end
%% Only use rats with 100 samples
[catStim,catBase] = deal(cell(1,2));
for ii = 1:size(thisStim,1)
    for jj = 1:size(thisStim,2)
        if size(thisStim{ii,jj},1)>=100 && ~any(any(isnan(thisStim{ii,jj})))
            this = thisStim{ii,jj}(randperm(size(thisStim{ii,jj},1),100),:);
            catStim{ii} = [catStim{ii};this]; 
            this = thisBase{ii,jj}(randperm(size(thisBase{ii,jj},1),100),:);
            catBase{ii} = [catBase{ii};this];
        end
    end
end
%% Build models using LOO - uses the fact that each animal has 100 samples
for ii = 1:2
    for jj = 1:size(catStim{ii},1)/100
        % Get indices
        testInds = (jj-1)*100+1:jj*100;
        trainInds = 1:size(catStim{ii},1);
        trainInds(testInds) = [];
        % Grab data
        testX = [catStim{ii}(testInds,:);catBase{ii}(testInds,:)];
        testY = [ones(100,1);zeros(100,1)];
        trainX = [catStim{ii}(trainInds,:);catBase{ii}(trainInds,:)];
        trainY = [ones(((size(catStim{ii},1)/100)-1)*100,1);...
            zeros(((size(catStim{ii},1)/100)-1)*100,1)];
        % Build and test model
        cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',[]);
        [~,~,~,~,acc{ii,jj},hist{ii,jj}] = lassoNet(trainX,trainY,...
            'binomial','deviance',1,10,1,cfg);
        % Run single features
        for k = 1:216
            mdl = fitglm(trainX(:,k),trainY,'distribution','binomial');
            pred = predict(mdl,testX(:,k));
            [~,~,~,a{ii,jj}(k)] = perfcurve(testY,pred,1);
            coeff{ii,jj}(k) = table2array(mdl.Coefficients(2,1));
        end
    end
end
save(['~/data/baseStim',num2str(n),'.mat'],'acc','hist','a','coeff')