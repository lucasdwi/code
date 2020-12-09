%% Compare coherence and power
% Calculate % change in features from baseline
for ii = 1:5
    percentChangeLSD{ii} =  (lsd{ii,2}-mean(lsd{ii,1}))./abs(mean(lsd{ii,1}));
    percentChangeSal{ii} = (saline{ii,2}-mean(saline{ii,1}))./abs(mean(saline{ii,1}));
end
% meanInds = [2:4,6:7];
% for ii = 1:5
%     percentChangeLSD{ii} = thisLSDDiff{ii}./abs(mean(lsdData{1}{meanInds(ii),1},1));
% %     percentChangeSal{ii} = thisBaseDiff{ii}./mean(allData{1}(meanInds(ii),1),1);
% end
deltaLSD = cat(1,percentChangeLSD{:});
deltaSal = cat(1,percentChangeSal{:});
% Set up logic to determine given a coherence feature index, what sites and
% freqeuency is involved
sites = {'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'};
freqs = {'d','t','a','b','lg','hg'};
cmbs = nchoosek(1:8,2);
powerInds = reshape(1:48,6,8)';
cohInds = lsdSortInd(lsdSortInd>=49);
figure
hold on
n = 10;
for ii = 1:n
    cohInd = cohInds(ii);
    this = cohInd-48;
    thisF = rem(this,6);
    if thisF == 0
        f = 6;
    else
        f = thisF;
    end
    thisCmb = ceil(this/6);
    s1 = powerInds(cmbs(thisCmb,1),f);
    s2 = powerInds(cmbs(thisCmb,2),f);
    % Plot both power features
    scatter([ii ii],[mean(deltaLSD(:,s1)) mean(deltaLSD(:,s2))].*100,'r')
    scatter(ii,mean(deltaLSD(:,cohInd)).*100,'k')
end
plot([0 n],[0 0],'--k')
set(gca,'xtick',1:n,'xticklabel',{feat{cohInds(1:n)}})
xtickangle(45)
legend({'Power','Coherence'})
ylabel('% change in stim from baseline')
title('LSD stim top coherence features')