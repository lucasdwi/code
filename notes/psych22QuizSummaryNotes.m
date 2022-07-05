ids = unique(q1(:,1));
scores = zeros(numel(ids),max(q1(:,2)));
for ii = 1:numel(ids)
    scores(ii,1:sum(q1(:,1)==ids(ii))) = cat(1,q1(q1(:,1)==ids(ii),3));
    for k = 2:max(q1(:,2))
        if scores(ii,k) == 0
            scores(ii,k) = NaN;
        end
    end
end
figure
hold on
for ii = 1:numel(ids)
   plot(scores(ii,:)./max(q1(:,3)),'.-k')
end
xlim([0.5 max(q1(:,2))])
set(gca,'XTick',1:1:max(q1(:,2)))
title(num2str(sum(~isnan(scores(:,2)))))
%%
for ii = 1:numel(ids)
    for k = 2:max(q1.Var2)
        if isnan(scores(ii,k))
            scores(ii,k) = scores(ii,k-1);
        end
    end
end
mean(scores)
%%
figure
hold on
violin(scores./18)
plotSpread(scores./18)