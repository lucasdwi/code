%% Pull out NACL theta
naclInd = find(T.Base.shellStrict==1);
naclNonInd = find(T.Base.shellStrict==0);
for ii = 1:length(naclInd);
    naclStrict(ii,:) = table2array(T.Base(naclInd(ii),18));
end
for ii = 1:length(naclNonInd);
    naclNon(ii,:) = table2array(T.Base(naclNonInd(ii),18));
end
naclAvg = mean(naclStrict); naclStd = std(naclStrict);
naclNonAvg = mean(naclNon); naclNonStd = std(naclNon);
%% Pull out 
lgamInd = find(T.Base.coreRespond==1);
lgamNonInd = find(T.Base.coreRespond==0);
for ii = 1:length(lgamInd);
    lgamResp(ii,:) = table2array(T.Base(lgamInd(ii),46));
end
for ii = 1:length(lgamNonInd);
    lgamNon(ii,:) = table2array(T.Base(lgamNonInd(ii),46));
end
lgamAvg = mean(lgamResp); lgamStd = std(lgamResp);
lgamNonAvg = mean(lgamNon); lgamNonStd = std(lgamNon);