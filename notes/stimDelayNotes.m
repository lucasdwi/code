for ii = 1:2
    % Find stop indices and times
    stopInd = logicFind(1,round(diff(eventTs.t{ii}),3)>0.1,'==');
    stopTime{ii} = eventTs.t{ii}(stopInd);
    % Get starts by going to the next index past stops
    startTime{ii} = eventTs.t{ii}(stopInd+1);
    % Add in last time
    stopTime{ii} = [stopTime{ii};eventTs.t{ii}(end)];
    % Add in first start
    startTime{ii} = [eventTs.t{ii}(1);startTime{ii}];
end
%%
startTimeDiff = startTime{2}-startTime{1};
stopTimeDiff = stopTime{2}-stopTime{1};

figure
plot(startTimeDiff,'.b')

figure
plot(stopTimeDiff,'.r')
%% Get means
inds = [1,10;11,510;511,1010;1011,1510;1511,2010;2011,2510;2511,3010];
for ii = 1:size(inds,1)
   mStart(ii) = mean(startTimeDiff(inds(ii,1):inds(ii,2)));
   mStop(ii) = mean(stopTimeDiff(inds(ii,1):inds(ii,2)));
end
