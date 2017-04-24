function [x,y,AUC] = cycleTest(fits,lambda,testX,testY)
%%
for ti = 1:size(testX,1)
    disp(num2str(ti))
    for r = 1:size(testY{ti,1},2)
    % Cycle test sets
    [predY] = cvglmnetPredict(fits{lambda{1,1}.bestLambdaInds(r),r},testX{ti},'lambda_1se','response');
    [x{ti}(:,r),y{ti}(:,r),~,AUC(r,ti)] = perfcurve(testY{ti}(:,r),predY,1);
    end
end