eventTimes 
%%
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\';
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\fNames');
%%
rest = cell(12,4); binge = cell(12,4);
cd(sdir)
for ci = 1:length(fNames)
    for ri = 1:size(fNames{ci},1)
    load([fNames{ci}{ri},'.mat'],'eventTs')
    % Get start and stop, and total time for rest and binge
    rest{ri,ci} = [eventTs.t{logicFind('Rest (Start)',eventTs.label,'==')},eventTs.t{logicFind('Rest (End)',eventTs.label,'==')}];
    rest{ri,ci}(:,3) = rest{ri,ci}(:,2)-rest{ri,ci}(:,1);
    binge{ri,ci} = [eventTs.t{logicFind('Binge (Start)',eventTs.label,'==')},eventTs.t{logicFind('Binge (End)',eventTs.label,'==')}];
    binge{ri,ci}(:,3) = binge{ri,ci}(:,2)-binge{ri,ci}(:,1);
    % Get toal time for rest and binge
    totRest(ri,ci) = sum(rest{ri,ci}(:,3));
    totBinge(ri,ci) = sum(binge{ri,ci}(:,3));
    % Get total time for recording
    total(ri,ci) = eventTs.t{1,1}(end); 
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
anovaBox(firstBingeNorm,group,'Binge Latency','% of trial')

%% Smooth mastertimeline
% Percent of animals bingeing
bingeing = sq(nanmean(masterTimeline,1));
figure
plot(bingeing)
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
   plot(masterTime,filtTime(:,ii))
   subplot(1,2,2)
   hold on
   plot(masterTime,bingeing(:,ii))
end
%% Compare binge size to time spent binging
%[r,p] = corrcoef(bingeSizes(~isnan(bingeSizes(:,3)),3),totBinge(~isnan(totBinge(:,3)),3))
figure
for ii = 1:4
   hold on
   plot(bingeSizes(:,ii),totBinge(:,ii),'.','markersize',10)
end
xlabel('Binge Size (gms)'); ylabel('Time Bingeing (s)')
%% Average speed
bingeSpeed = bingeSizes./totBinge;
% Run ANOVA
[p,tbl,stats] = anova1(bingeSizes,group,'off');
% Run multiple corrections (Bonferroni) on all pairs of means
[mctbl] = multcompare(stats,'CType','bonferroni','display','off');
% Combine group labels for sigstar
for ii = 1:size(mctbl,1)
    sigGroups{ii} = [mctbl(ii,1),mctbl(ii,2)];
end
% Set-up boxplot
figure; 
boxplot(bingeSizes,group)
hold on
% Plot means as red squares
for ii = 1:length(group)
   plot(ii,stats.means(ii),'rs') 
end
title('Binge Speed across Conditions')
ylabel('Binge Speed (gm/sec)');
% Find indices of mean comparisons with significant differences
sigInds = logicFind(0.05,mctbl(:,6),'<=');
% Plot significance bars with sigstar
sigstar(sigGroups(sigInds),mctbl(sigInds,6))
%% Ephys T-Tests