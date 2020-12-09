%% Frequency networks
% Raw inputs
freq = {'d','t','a','b','lg','hg'};
sites = {'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC','rNAcC'};
chan = numel(sites);
% Set input data to use as edges (coh) and node (power) sizes
thisData = aLSDsingle'.*100;
% Set random data to use to set significance threshold - same second
% dimension
thisDataR = aLSDsingleP'.*100;
% What to use to indicate direction (sign)
thisSign = sign(mean(betaLSD,2));
% Get number of power features
powerN = chan*numel(freq);
%%
% Set vector of mean data
thisVectM = mean(thisData,1);
% Get combinations of channels
cmbs = nchoosek(1:chan,2);
% Run ttest2
[~,p] = ttest2(thisData,thisDataR);
% Simple bonferroni correction (slightly over-conseravative, higher false negative)
pAdj = p.*numel(thisVectM);
% Find scaling factor based on minimum significant for node size
% s = 1-min(abs(thisVectM(pAdj<=0.05)));
s = -50.9384; % manually set to match across lsd and saline data sets
% Set flexible scales for node and edge weights (so they min and max in visible scale)
pow = thisVectM(1:powerN);
coh = thisVectM(powerN+1:end);
powV = pow(pAdj(1:powerN)<=0.05);
cohV = coh(pAdj(powerN+1:end)<=0.05);
maxPow = max(powV); minPow = min(powV);
maxCoh = max(cohV); minCoh = min(cohV);


scaled = linspace(25,75,40);
test = scaled(realSize-25);
% figure
for jj = 1%:6
    adjMat = [];
    edgeSign = [];
    count = powerN+jj;
    c = 1;
    for ii = 1:size(cmbs,1)
        % Make sure edge is significant
        if pAdj(count) <= 0.05
            adjMat(cmbs(ii,1),cmbs(ii,2)) =  thisVectM(count);
            edgeSign(c) = thisSign(count);
            c=c+1;
        else
            adjMat(cmbs(ii,1),cmbs(ii,2)) = 0;
        end
        count = count+numel(freq);
    end
    adjMat = [adjMat;zeros(1,8)];
    % Compute graph with scale factor
    g = graph(adjMat,sites,'upper');
    % Set up node coordinates based on regular octagon
    x = [1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2,0,sqrt(2)/2];
    y = [0,sqrt(2)/2,1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2];
    % Set edge color based on edge sign (+ = red; - = blue)
    edgeColor = zeros(numel(edgeSign),3);
    edgeColor(logicFind(1,edgeSign,'=='),:) = repmat([1 0 0],sum(edgeSign==1),1);
    edgeColor(logicFind(-1,edgeSign,'=='),:) = repmat([0 0 1],sum(edgeSign==-1),1);
    % Set node color based on node sign (+ = red; - = blue)
    nodeSign = thisSign(jj:6:powerN)>=0;
    nodeColor = zeros(numel(nodeSign),3);
    nodeColor(logicFind(1,nodeSign,'=='),:) = repmat([1 0 0],sum(nodeSign==1),1);
    nodeColor(logicFind(0,nodeSign,'=='),:) = repmat([0 0 1],sum(nodeSign==0),1);
    % Check node significance
    sigNode = pAdj(jj:6:powerN) <= 0.05;
    % If not significant, set to white and size 10
    nodeColor(~sigNode,:) = repmat([1 1 1],sum(~sigNode),1);
    theseNodes = abs(thisVectM(jj:6:48));
    realSize = round(theseNodes);
    theseNodes(sigNode) = scaled(realSize(sigNode)-50);
    theseNodes(~sigNode) = 25;
    % Plot
    figure('position',[895,400,667,571])
% %     subplot(2,3,jj)
    h = plot(g,'XData',x,'YData',y,'markersize',theseNodes,'linewidth',abs(g.Edges.Weight)+s,'edgecolor',edgeColor,'nodecolor',nodeColor);
    title(freq{jj})
    axis off
    axis square
    print(['LSDvSal_stimImag_logLoo_all-216',freq{jj},'.eps'],'-dwinc','-painters')
end