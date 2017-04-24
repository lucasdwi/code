%% Ephys Enet
files = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper2\percBingeRand\'},{'.mat'});
%files = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\subsetLasso\'},{'.mat'});
%files = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\ephys4foldEnet\'},{'enet4fold'});
%files = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\behavior3foldEnet\'},{'enet3foldBehavior'});
beta = [];
masterErr = []; bestModel = []; masterRawErr = [];
for fi = 1:size(files{1},1)
    load(files{1}(fi).name,'allBeta','allLambda','allAlpha')%,'cvFitsArray')
    for s = 1:size(allBeta,2)
       beta = cat(1,beta,allBeta{s}.betas);
    end
    for s = 1:size(allLambda,2)
        masterErr = [masterErr;allLambda{s}.allErr];
%         bestModel = cat(3,bestModel,[allAlpha{s}.bestAlpha;allLambda{s}.bestLambda;allLambda{s}.minLamErr]);
    end
%     for a = 1:size(cvFitsArray,2)
%         for r = 1:size(cvFitsArray{a},1)
%             masterRawErr = [masterRawErr;cvFitsArray{a}{r,1}.cvm];
% %             masterRawErr = [masterRawErr;[cvFitsArray{a}{r,1}.cvm,cvFitsArray{a}{r,2}.cvm,cvFitsArray{a}{r,3}.cvm,cvFitsArray{a}{r,4}.cvm]];
%         end
%     end
end
%%
bestModelMean = mean(bestModel,3);
bestModelSd = std(bestModel,[],3);
%%
masterErrAvg = mean(masterErr,1);
masterErrSd = std(masterErr,[],1);
%% Conduct t-test
% Subtract 0.5 from all data to center at zero
centerMaster = masterErr - 0.5;
%%
rawErr = [];
for c = 1:size(cvFitsArray{1},2)
    for r = 1:size(cvFitsArray{1},1)
        rawErr = [rawErr;[cvFitsArray{1}{r,1}.cvm,cvFitsArray{1}{r,2}.cvm,cvFitsArray{1}{r,3}.cvm,cvFitsArray{1}{r,4}.cvm]];
    end
end
%% One structure 
masterErr = [];
for ii = 1:size(allLambda,2)
    masterErr = [masterErr;normLambda{1,ii}.allErr];
end

    