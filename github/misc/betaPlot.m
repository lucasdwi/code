%% Plots beta values
ys = T.Base.Properties.VariableNames(8:57);
for iy = 1:length(ys)
    ys{iy} = strrep(ys{iy},'_',' ');
end
%%
for b = 1:size(allDev,2)
    % Find indices of good betas
    goodI = find(~isnan(mBeta(b,:)));
    figure
    scatter((mBeta(b,goodI)),1:length(goodI),'.k')
    set(gca,'YTick',1:length(goodI),'YTickLabel',ys(goodI));
    hold on;
    % Plot error bars
    for h = 1:length(goodI)
        plot([mBeta(b,goodI(h))-stdBeta(b,goodI(h)),mBeta(b,goodI(h))+stdBeta(b,goodI(h))],[h,h],'k');
    end
    % Plot line at x = 0
    plot([0,0],get(gca,'YLim'),'color',[0.3 0.3 0.3]);
end
%% One plot
% Get indices with at least one non-NaN in column
inds = find(~isnan(nanmean(mBeta,1)));
% Get titles for plots
ts = T.Base.Properties.VariableNames(2:5);
figure; hold on;
for sp = 1:4
    subplot(2,2,sp); %subaxis(2,2,sp,'SpacingHoriz',0.3,'SpacingVert',0.1,'Padding',0,'Margin',0.1);
    goodI = find(~isnan(mBeta(sp,:)));
    plot(mBeta(sp,goodI),goodI,'.k');
    hold on;
    for h = 1:length(find(~isnan(mBeta(sp,:))));
        plot([mBeta(sp,goodI(h))-stdBeta(sp,goodI(h)),mBeta(sp,goodI(h))+stdBeta(sp,goodI(h))],[goodI(h),goodI(h)],'k');
    end
    ylim([0 50]); xlim([-2 3]);
    plot([0 0],get(gca,'YLim'),'color',[.8 .8 .8]);
    set(gca,'YTick',goodI,'YTickLabel',ys(goodI));
    title(ts{sp})
end
%% Multiple model iteration beta heat-map
% Get titles
ts = T.Base.Properties.VariableNames(2:5);
for m = 1:4
    thisOR = [];
    for ii = 1:length(masterOR)
        % Concatenate all iterations of model
        thisOR = horzcat(thisOR,masterOR{1,ii}(m,:)');
    end
    goodI = find(~isnan(nanmean(thisOR,2)));
    figure;
    h = pcolor([thisOR nan(size(thisOR,1),1); nan(1,size(thisOR,2)+1)]);
    % Set x and y ticks and labels
    xlabel('Iteration');
    set(gca,'XTick',[1.5,2.5,3.5,4.5,5.5],'XTickLabel',[1,2,3,4,5],'YTick',goodI+.75,'YTickLabel',ys(goodI));
    set(gca,'FontSize',8);
    title(ts(m))
end
%% Plot beta and p(t)
% Extract and threshold p values from T-Tests
pMat = cell2mat(p);
pT(:,1:2) = pMat(101:150,1:2);
pT(:,3:4) = pMat(201:250,1:2);
pT(pT>=0.05) = NaN;
% Combine p values and ORs
for ii = 1:4
    POR(:,ii*2-1:ii*2) = horzcat(mOR(ii,:)',pT(:,ii));
end
% Get indices that have either p or OR values
goodI = find(~isnan(nanmean(POR,2)));
goodPOR = POR(goodI,:);
% NaN pad
goodPOR = horzcat(goodPOR,NaN(size(goodPOR,1),1));
goodPOR = vertcat(goodPOR,NaN(1,size(goodPOR,2)));
%% Sublot mOR and p values
oIs = [1:2:9]; pIs = [2:2:9];
goodOR = goodPOR; goodOR(:,pIs) = NaN;
goodP = goodPOR; goodP(:,oIs) = NaN;
figure; 
subplot(1,2,1)
im1 = pcolor(goodOR);
set(im1.Parent,'XTick',[2,4,6,8],'XTickLabel',{'Shell Respond','Shell Strict','Core Respond','Core Strict'},'XTickLabelRotation',45,'YTick',[8,24],'YTickLabel',{'Power','Coherence'},'YTickLabelRotation',90)
hold on;
set(im1,'alphadata',~isnan(goodOR));
% Horizontal line splitting power from coherence variables
plot([0 9],[13 13],'k','LineWidth',2);
% Vertical lines splitting each model
plot([3 3],[0 33],'k','LineWidth',1); plot([5 5],[0 33],'k','LineWidth',1); plot([7 7],[0 33],'k','LineWidth',1);

subplot(1,2,2)
im2 = pcolor(goodP);
set(im2,'alphadata',~isnan(goodP));
set(im2.Parent,'XTick',[2,4,6,8],'XTickLabel',{'Shell Respond','Shell Strict','Core Respond','Core Strict'},'XTickLabelRotation',45,'YTick',[8,24],'YTickLabel',{'Power','Coherence'},'YTickLabelRotation',90)
hold on;
set(im1,'alphadata',~isnan(goodOR));
% Horizontal line splitting power from coherence variables
plot([0 9],[13 13],'k','LineWidth',2);
% Vertical lines splitting each model
plot([3 3],[0 33],'k','LineWidth',1); plot([5 5],[0 33],'k','LineWidth',1); plot([7 7],[0 33],'k','LineWidth',1);

    