%%
sdir = 'D:\paper2\mat\';
[fNames(:,1)] = fileSearch(sdir,'Base','in');
[fNames(:,2)] = fileSearch(sdir,'Dep24','in');
[fNames(1:9,3)] = fileSearch(sdir,'Dep48','in');
[fNames(1:9,4)] = fileSearch(sdir,'Chow','in');
%%
rest = cell(12,4); binge = cell(12,4); start = cell(12,4);
% cd(sdir)
for ri = 1:size(fNames,1)
    for ci = 1:size(fNames,2)
        if ~isempty(fNames{ri,ci})
            load([fNames{ri,ci}],'eventTs')
            % Get absolute start time
            start{ri,ci} = eventTs.t{logicFind('food',eventTs.label,'==')};
            % Get all rest starts and stops
            rest{ri,ci} = [eventTs.t{logicFind('Rest (Start)',eventTs.label,'==')},eventTs.t{logicFind('Rest (End)',eventTs.label,'==')}];
            % Normalize
            rest{ri,ci} = rest{ri,ci}-start{ri,ci};
            % Get length per epoch
            rest{ri,ci}(:,3) = rest{ri,ci}(:,2)-rest{ri,ci}(:,1);
            % Get all binge starts and stops
            binge{ri,ci} = [eventTs.t{logicFind('Binge (Start)',eventTs.label,'==')},eventTs.t{logicFind('Binge (End)',eventTs.label,'==')}];
            % Normalize
            binge{ri,ci} = binge{ri,ci}-start{ri,ci};
            % Get length per epoch
            binge{ri,ci}(:,3) = binge{ri,ci}(:,2)-binge{ri,ci}(:,1);
            % Get toal time for rest and binge
            totRest(ri,ci) = sum(rest{ri,ci}(:,3));
            totBinge(ri,ci) = sum(binge{ri,ci}(:,3));
            % Get total time for recording
            total(ri,ci) = eventTs.t{1,1}(end);
        end
    end
end
% Get number of rests and binges
nRest = cellfun(@(sz) size(sz,1),rest);
nBinge = cellfun(@(sz) size(sz,1),binge);
% Set 0 to NaN
totRest(totRest == 0) = NaN;
totBinge(totBinge == 0) = NaN;
nRest(nRest == 0) = NaN;
nBinge(nBinge == 0) = NaN;
% Get average time per epoch
avgRest = totRest./nRest;
avgBinge = totBinge./nBinge;
%% Construct normalized binge timelines
timeline = cell(12,4);
% Use largest file as master normalized time steps
% Percent of length per half second
timeStep = 1/max(max(total));
masterTime = 0:timeStep:1;
masterTimeline = zeros(12,length(masterTime),4);
for ci = 1:size(total,2)
    for ri = 1:size(total,1)
        % Set up timeline
        timeline{ri,ci} = [0:0.5:total(ri,ci)];
        if size(timeline{ri,ci},2) > 1
            % Get normalized timeline
            timeline{ri,ci}(2,:) = timeline{ri,ci}(1,:)/total(ri,ci);
            % Insert 1s for binges
            for bi = 1:nBinge(ri,ci)
                % Each animals' timeline
                timeline{ri,ci}(3,nearest_idx3(binge{ri,ci}(bi,1),timeline{ri,ci}(1,:),-1):nearest_idx3(binge{ri,ci}(bi,2),timeline{ri,ci}(1,:))) = 1;
                % Master timeline with same normalized time axis
                masterTimeline(ri,nearest_idx3(binge{ri,ci}(bi,1)/total(ri,ci),masterTime(1,:)):nearest_idx3(binge{ri,ci}(bi,2)/total(ri,ci),masterTime(1,:)),ci) = 1;
            end
            % Get normalized time to first binge
            firstBingeNorm(ri,ci) = binge{ri,ci}(1,1)/total(ri,ci);
        else
            % Set to NaN if empty
            masterTimeline(ri,:,ci) = NaN;
        end
    end
end
%% Run ANOVA on time to first binge
% Get rid of zeros
firstBingeNorm(firstBingeNorm == 0) = NaN;
% Run ANOVA
x = reshape(firstBingeNorm,1,48);
x(isnan(x)) = [];
group = [repmat(1,1,size(fNames{1},2)),repmat(2,1,size(fNames{2},2)),repmat(3,1,size(fNames{3},2)),repmat(4,1,size(fNames{4},2))];
anovaBox(x,group,'Binge Latency','% of trial')
%%
figure
plot([repmat(1,12,1),repmat(2,12,1),repmat(3,12,1),repmat(4,12,1)],firstBingeNorm,'.k')
hold on
plot([1,2,3,4],[mean(firstBingeNorm(:,1)),mean(firstBingeNorm(:,2)),mean(firstBingeNorm(:,3),'omitnan'),mean(firstBingeNorm(:,4),'omitnan')],'rs')
%% Smooth mastertimeline
% Percent of animals bingeing
bingeing = squeeze(nanmean(masterTimeline,1));
figure
plot(bingeing)
% Split into ten bins (~10%)
inds = nearest_idx3([0.1:0.1:1],masterTime);
inds(:,2) = [1;875;1748;2620;3493;4366;5239;6112;6984;7857];
for ii = 1:length(inds)
    for j = 1:4
        m(ii,j) = mean(bingeing(inds(ii,2):inds(ii,1),j));
    end
end
%% Plot
figure
hold on
for ii = 1:4
    plot([0.1:0.1:1].*100,m(:,ii).*100)
end
xlabel('% of Session')
ylabel('Average % of Animals Bingeing')
title('% of Animals Bingeing across Conditions')
legend({'Base','Dep24','Dep48','Chow'},'location','northeast')
%% Find latency to 50% bingers
for ii = 1:4
    inds = logicFind(0.5,bingeing(:,ii),'>=');
    firsts(ii) = masterTime(inds(1));
end
%% Find latency to max % bingers
maxPerc = max(bingeing);
for ii = 1:4
    inds = logicFind(maxPerc(ii),bingeing(:,ii),'==');
    latMax(ii) = masterTime(inds(1));
end
%% 'Survival' Curve
figure
for ii = 1:4
    ecdf(firstBingeNorm(:,ii))
    hold on
end
legend({'base','dep24','dep48','chow'},'location','southeast')
title('First Binge')
xlabel('% of Session')
ylabel('Cumulative Density')
cmbs = nchoosek(1:4,2);
for ii = 1:size(cmbs,1)
    [h(ii),cmbs(ii,3)] = kstest2(firstBingeNorm(:,cmbs(ii,1)),firstBingeNorm(:,cmbs(ii,2)));
end
%%
sigma = 75;
sz = 8728;
x = linspace(-sz/2,sz/2,sz);
gFilt = exp(-x.^2/(2*sigma^2));
gFilt = gFilt/sum(gFilt);

figure
clear filtTime
for ii = 1:4
    filtTime(:,ii) = conv(gFilt,bingeing(:,ii),'same');
    subplot(1,2,1)
    hold on
    plot(masterTime,filtTime(:,ii));
    subplot(1,2,2)
    hold on
    plot(masterTime,bingeing(:,ii))
end
legend({'base','dep24','dep48','chow'},'location','northeast')
%% Compare binge size to time spent binging
% Set up aoc table
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
aTab = [bingeCal(:,1);bingeCal(:,2);bingeCal(:,3);bingeCal(:,4)];
aTab(:,2) = [totBinge(:,1);totBinge(:,2);totBinge(:,3);totBinge(:,4)];
aTab(isnan(aTab(:,1)),:) = [];
group = {'base','base','base','base','base','base','base','base','base','base','base','base','dep24','dep24','dep24','dep24','dep24','dep24','dep24','dep24','dep24','dep24','dep24','dep24','dep48','dep48','dep48','dep48','dep48','dep48','dep48','dep48','dep48','chow','chow','chow','chow','chow','chow','chow','chow',};
aoctool(aTab(:,1),aTab(:,2),group)
%[r,p] = corrcoef(bingeSizes(~isnan(bingeSizes(:,3)),3),totBinge(~isnan(totBinge(:,3)),3))
%%
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
figure
for ii = 1:4
    hold on
    h{ii} = plot(bingeCal(:,ii),bingeTime(:,ii),'.','markersize',10);
    % Plot regression lines
    lsline
end
xlabel('Binge Size (kCal)'); ylabel('Time Bingeing (s)')
legend([h{1},h{2},h{3},h{4}],{'Base','Dep24','Dep48','Chow'},'location','northwest')
title('Binge Size vs. Time Spent Bingeing')
%% Average speed
bingeSpeed = bingeCal./bingeTime;
% Run ANOVA
[p,tbl,stats] = anova1(bingeSpeed,[],'off');
% Run multiple corrections (Bonferroni) on all pairs of means
[mctbl] = multcompare(stats,'CType','bonferroni','display','off');
% Combine group labels for sigstar
for ii = 1:size(mctbl,1)
    sigGroups{ii} = [mctbl(ii,1),mctbl(ii,2)];
end
% Set-up boxplot
figure;
plot(repmat(1:4,12,1),bingeSpeed,'.k')
xlim([0 5])
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'Base','Dep24','Dep48','Chow'})
%boxplot(bingeSpeed,{'base','dep24','dep48','chow'})
hold on
% Plot means as red squares
for ii = 1:size(bingeSpeed,2)
    plot(ii,stats.means(ii),'rs')
end
title('Binge Speed across Conditions')
ylabel('Binge Speed (kCal/sec)');
% Find indices of mean comparisons with significant differences
sigInds = logicFind(0.05,mctbl(:,6),'<=');
% Plot significance bars with sigstar
sigstar(sigGroups(sigInds),mctbl(sigInds,6))
%% Binge size cals
[p,tbl,stats] = anova1(bingeCal,[],'off');
[mctbl] = multcompare(stats,'CType','bonferroni','display','off');
for ii = 1:size(mctbl,1)
    sigGroups{ii} = [mctbl(ii,1),mctbl(ii,2)];
end
figure;
plot(repmat(1:4,12,1),bingeCal,'.k')
xlim([0 5])
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'Base','Dep24','Dep48','Chow'})
hold on
for ii = 1:size(bingeCal,2)
    plot(ii,stats.means(ii),'rs')
end
title('Binge Size across Conditions')
ylabel('Calories (kCal)')
sigInds = logicFind(0.05,mctbl(:,6),'<=');
sigstar(sigGroups(sigInds),mctbl(sigInds,6))