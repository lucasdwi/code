function permTest(dist1,dist2,n)
%% Check if input distributions are columns, else transpose
if ~iscolumn(dist1)
    dist1 = dist1';
end
if ~iscolumn(dist2)
    dist2 = dist2';
end
%% Get observed difference/test statistic
dObv = distES(dist1,dist2);
% obv = mean(dist1)-mean(dist2);
% Get sample sizes of distributions
n1 = numel(dist1);
n2 = numel(dist2);
% Pool distributions
pool = [dist1;dist2];
nPool = numel(pool);
% Cycle through n iterations of resampling and calculating test-statistic
for permIter = 1:n
    perm1 = pool(randperm(nPool,n1));
    perm2 = pool(randperm(nPool,n2));
%     d(permIter) = mean(perm1) - mean(perm2);
    d(permIter) = distES(perm1,perm2);
end
%%
(sum(abs(d)>=dObv)+1)/n;

