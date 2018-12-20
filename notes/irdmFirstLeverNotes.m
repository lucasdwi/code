%% Get stim times and separate files (manual)
oldEventTs = eventTs;
oldLFPTs = LFPTs;
%% Grab stims and data
eventTs = oldEventTs;
LFPTs = oldLFPTs;
% Find stims
stimChan = 15;

% Separate stim and clean data
stimInd = []; stim = []; stimStop = []; stimStart = [];
% stimInd = logicFind(1,~cellfun(@isempty,eventTs.t),'==');
stim = eventTs.t{1,stimChan};
% Get all stops
stimStop = stim(logicFind(0,round(diff(stim)),'~='));
% Add last stop since above algorithm skips
stimStop = [stimStop;stim(end)];
% Go 0.5 seconds before each stop
stimStop = stimStop;
% Get starts using the next index after stops from stim
stimStart = stim(logicFind(0,round(diff(stim)),'~=')+1);
% Add first start since above algorithm skips
stimStart = [stim(1);stimStart];
% Use stim start and stops to determine data start and stops; also include
% first/last data point - assumes that no stim starts or ends a recording
dataStart = [LFPTs.tvec(1);LFPTs.tvec(nearest_idx3(stimStop,LFPTs.tvec)+1)'];
dataStop = [LFPTs.tvec(nearest_idx3(stimStart,LFPTs.tvec)-1)';LFPTs.tvec(end)];

eventTs.label{1} = 'Inter (Start)';
eventTs.label{2} = 'Inter (End)';
eventTs.t{1} = dataStart;
eventTs.t{2} = dataStop;
% Grab data
LFPTs.data = LFPTs.data(17:24,:);
LFPTs.label = LFPTs.label(17:24);
eventTs.label = eventTs.label([1:2,49:54]);
eventTs.t = eventTs.t([1:2,49:54]);
%% Set up inputs
nFilt = [57 63];
dsf = 5;
thresh = 0.5;
onset = 0.25;
offset = 0.25;
foi = [1 1 100];
len = [0.25,0.5,0.75,1:.5:5];
fNames = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\irdm\'...
    'stimMat\'],'.mat');
stop = 0;
for ii = 1:size(fNames,2)
    disp(num2str(ii))
    for jj = 1:size(len,2)
        load(fNames{ii})
        LFPTs.data = LFPTs.data(2:8,:);
        [LFPTs.data] = filter60(LFPTs,adfreq,0);
        [LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq);
        eoi = {'inter',[0 len(jj)]};
        minInt = max(diff(cell2mat(eoi(:,2)),1,2));
        [LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,...
            offset,minInt,adfreq);
        [clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
        if ~isempty(trls{1})
            nTrial(ii,jj) = size(trls{1}.trial,3);
        else
            nTrial(ii,jj) = 0;
        end
    end
end
%%
figure
for jj = 1:17
    hold on
    for ii = 1:7
        plot(trls{1,1}.time{jj},trls{1,1}.trial(ii,:,jj))
    end
    title(num2str(jj))
    pause(2.5)
    cla
end
%%
fNames = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\irdm\mat\','.mat');
for fi = 1:size(fNames,2)
    disp([num2str(fi),' of ',num2str(size(fNames,2))])
    load(fNames{fi})
    % Create event 'trial' that corresponds to every time both lights turn on
    % at the same time; there should be 10x5 instances
    ind = zeros(size(eventTs.t{1,5},1),1);
    for ii = 1:size(eventTs.t{1,5},1)
        this = logicFind(eventTs.t{1,5}(ii),eventTs.t{1,6},'==');
        if ~isempty(this)
            ind(ii) =  this;
        end
    end
    % Check number of trials
    if sum(ind~=0) ~= 50
        error(['Only ',num2str(sum(ind~=0)),' trials were found.'])
    end
    % Create event label and timestamps
    eventTs.label{1,end+1} = 'trial';
    eventTs.t{1,end+1} = eventTs.t{1,5}(ind~=0);
    % For each trial, determine which lever was pressed first
    lever1 = eventTs.t{1,1}(nearest_idx3(eventTs.t{1,7},eventTs.t{1,1},1));
    lever2 = eventTs.t{1,2}(nearest_idx3(eventTs.t{1,7},eventTs.t{1,2},1));
    % Set-up difference vectors
    diffvec = [lever1-eventTs.t{1,7},lever2-eventTs.t{1,7}];
    % NaN negative values
    diffvec(diffvec<0) = NaN;
    % Find minumum time difference
    [~,leverInd] = min(diffvec,[],2,'omitnan');
    % Create variable for which lever was pressed first in each trial; to be
    % used as a dependent variable in prediction
    firstLever = eventTs.label(leverInd);
    save(fNames{fi},'LFPTs','adfreq','eventTs','pl2','firstLever')
end
%% Collate data and get used dependent variable (firstLever)
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\irdm\processed\'],{'.mat'},{'pow','coh'},'trl','rel');
for ii = 1:size(files{1},2)
    load(['C:\Users\Pythia\Documents\GreenLab\data\irdm\mat\',...
        files{1}{ii}(1:end-10),'.mat'])
    % Find the trials that match the samps used; accounting for downsample
    % factor of 5
    inds = nearest_idx3(LFPTs.tvec(samp{1}{ii}(:,2)*5),...
        eventTs.t{1,7});
    lever{ii} = firstLever(inds);
end
%% Separate data by lever pressed
catData = cell(1,2);
for ii = 1:size(lever,2)
    % Replace ILP with 1 and DLP with 0
    binLever{ii} = strcmp(lever{1,ii},'ILP');
    catData{1,1} = [catData{1,1};data{1,1}{ii}(binLever{ii},:)];
    catData{1,2} = [catData{1,2};data{1,1}{ii}(~binLever{ii},:)];
end
%% Build models
for ii = 1:100
    % Take out 15 of each for testing
    testInds1 = randperm(66,15);
    testInds2 = randperm(68,15);
    testX{ii} = [catData{1,1}(testInds1,:);catData{1,2}(testInds2,:)];
    
    trainInds1 = ~ismember(1:66,testInds1);
    trainInds2 = ~ismember(1:68,testInds2);
    % Impute up
    newX1 = normOversample(catData{1,1}(trainInds1,:),100);
    newX2 = normOversample(catData{1,2}(trainInds2,:),100);
    trainX{ii} = [catData{1,1}(trainInds1,:);newX1;catData{1,2}(trainInds2,:);newX2];
    %     trainX{ii} = [catData{1,1}(trainInds1,:);catData{1,2}(trainInds2,:)];
    %     mdl{ii} = fitglm(trainX{ii},[ones(66-15,1);zeros(68-15,1)],...
    %         'distribution','binomial');
    mdl{ii} = fitglm(trainX{ii},[ones(66-15+100,1);zeros(68-15+100,1)],...
        'distribution','binomial');
    prob = predict(mdl{ii},testX{ii});
    [x(ii,:),y(ii,:),~,a(ii)] = perfcurve([ones(15,1);zeros(15,1)],prob,...
        1,'TVals',0:1/99:1,'UseNearest',0);
    % Permuted
    prob = predict(mdl{ii},testX{ii}(randperm(30,30),:));
    [xRnd(ii,:),yRnd(ii,:),~,aRnd(ii,:)] = perfcurve([ones(15,1);...
        zeros(15,1)],prob,1,'TVals',0:1/99:1,'UseNearest',0);
end
%% Exmaple data
figure
hold on
plot(LFPTs.tvec,LFPTs.data(2,:))
for ii = 1:7
    plot(eventTs.t{1,ii},1+0.5*ii,'.')
    if ii == 7
        plot(eventTs.t{1,ii},1+0.5*ii,'ob')
    end
end
%%
figure
hold on
for ii = 1:8
    plot(LFPTs.tvec,LFPTs.data(ii,:)+ii-1)
end
% shapes = {'.r','.b','sqc','*y','ob','or','om'};
% for ii = 1:7
%     plot(eventTs.t{1,ii},-1,'.')
% end
%% Compare base to ddt stim
[ddt,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\irdm\processedStim\',{'DDT'},{'pow','coh'},'avg','rel');
[base,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\irdm\processedStim\',{'irdm'},{'pow','coh'},'trl','rel');
inds = {(5:8),(9:12),(13:14)};
for ii = 1:size(ddt{1},1)
    z(ii,:) = (ddt{1,1}{ii}-mean(cat(1,base{1}{inds{ii}}),1))./...
        std(cat(1,base{1}{inds{ii}}),[],1);
end
figure
h = pcolor(padarray(z,1,'post')');
set(h,'edgecolor','none');
colormap viridis
set(gca,'xtick',1.5:3.5,'xticklabel',[2,5,6])
xlabel('animal')
ylabel('feature')
title('Base vs. Stim (z-score)')
c = colorbar;
c.Label.String = 'z-score';