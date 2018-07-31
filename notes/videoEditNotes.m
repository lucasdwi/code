tic
% Create object to get number of frames
vIn = VideoReader('C:\Users\Pythia\Documents\GreenLab\data\paper2\I3BaseNov11.AVI');
nF = vIn.NumberOfFrames;
% Then recreate object to loop through
vIn = VideoReader('C:\Users\Pythia\Documents\GreenLab\data\paper2\I3BaseNov11.AVI');
% Create output object
vOut = VideoWriter('test.avi');
% Load trls structure for samples
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge\I3Base_2015-11-11_binge_vs_notbinge.mat', 'trls','LFPTs')
allSamp = 1:length(LFPTs.tvec);
feed = [];
for ii = 1:size(trls{1,1}.sampleinfo,1)
    feed = cat(2,feed,trls{1,1}.sampleinfo(ii,1):trls{1,1}.sampleinfo(ii,2));
end
notFeed = [];
for ii = 1:size(trls{1,2}.sampleinfo,1)
    notFeed = cat(2,notFeed,trls{1,2}.sampleinfo(ii,1):trls{1,2}.sampleinfo(ii,2));
end
noise = allSamp(~ismember(allSamp,[feed,notFeed]));
% Create event timeline, 0s = noise, 1s = feed, 2s = notFeed
event = zeros(1,length(allSamp));
event(feed) = 1;
event(notFeed) = 2;
% Set up colors
col = [0,255,0;255,0,0;255,255,255];
% Open output object
open(vOut)
% Set frame counter, k
k = 1;
% Loop trough all frames of vIn, and add element to corner
while k <= nF
    frame = readFrame(vIn);
    % Use previous stopInd + 1 as new startInd, unless k == 1, then
    % startInd = 1
    if k == 1
        startInd = 1;
    else
        startInd = stopInd + 1;
    end
    % Get nearest time stamp of end of frame
    stopInd = nearest_idx3(k*1/vIn.FrameRate,LFPTs.tvec);
    % Use majority event during frame for color
    dummy = event(startInd:stopInd);
    % If tie, deaults to first max
    [~,beh] = max([sum(dummy==0),sum(dummy==1),sum(dummy==2)]);
    % Color upper right square according to behavior
    frame(1:20,621:640,:) = cat(3,ones(20,20).*col(beh,1),ones(20,20).*col(beh,2),ones(20,20).*col(beh,3));
    writeVideo(vOut,frame);
    k = k +1;
end
close(vOut)
toc
%%
load('baseline500Each6000All50-50.mat')
animals = 1:12;
inds = 1:60;
inds = inds(~ismember(inds,pInds));
trainX = cat(1,each.trainX{1,animals(~ismember(animals,3))});
trainX = trainX(:,inds);
trainY = cat(1,each.trainY{1,animals(~ismember(animals,3))});
mdl = fitglm(trainX,trainY,'distribution','binomial');
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\I3Base_video_predict.mat')
prob = predict(mdl,x(:,inds));
[rocX,rocY,~,a] = perfcurve(y,prob,1);
%%
tic
% Create object to get number of frames
vIn = VideoReader('C:\Users\Pythia\Documents\GreenLab\data\paper2\I3BaseNov11.AVI');
nF = vIn.NumberOfFrames;
% Then recreate object to loop through
vIn = VideoReader('C:\Users\Pythia\Documents\GreenLab\data\paper2\I3BaseNov11.AVI');
% Create output object
vOut = VideoWriter('C:\Users\Pythia\Documents\GreenLab\data\paper2\I3BaseNove11Pred.AVI');
% Load trls structure for samples
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge\I3Base_2015-11-11_binge_vs_notbinge.mat', 'trls','LFPTs')

allSamp = 1:length(LFPTs.tvec);
pred = round(prob);
pred(pred==0) = 2;
feed = [];
c = 1;
for ii = 1:size(trls{1,1}.sampleinfo,1)
    feed(1,2000*ii-1999:ii*2000) = trls{1,1}.sampleinfo(ii,1):trls{1,1}.sampleinfo(ii,2);
    feed(2,2000*ii-1999:ii*2000) = pred(c);
    c = c+1;
end
notFeed = [];
c = 264;
for ii = 1:size(trls{1,2}.sampleinfo,1)
    notFeed(1,2000*ii-1999:ii*2000) = trls{1,2}.sampleinfo(ii,1):trls{1,2}.sampleinfo(ii,2);
    notFeed(2,2000*ii-1999:ii*2000) = pred(c);
    c = c+1;
end

event = zeros(1,length(allSamp));
event(feed) = 1;
event(notFeed) = 2;

eventPred = zeros(1,length(allSamp));
eventPred(feed(1,:)) = feed(2,:);
eventPred(notFeed(1,:)) = notFeed(2,:);

k = 1;
% Set up colors
col = [255,0,0;255,255,255;0,0,0];
% Open output object
open(vOut)
while k <= nF
    frame = readFrame(vIn);
    % Use previous stopInd + 1 as new startInd, unless k == 1, then
    % startInd = 1
    if k == 1
        startInd = 1;
    else
        startInd = stopInd + 1;
    end
    % Get nearest time stamp of end of frame
    stopInd = nearest_idx3(k*1/vIn.FrameRate,LFPTs.tvec);
    % Use majority event during frame for color
%     dummy = event(startInd:stopInd);
%     % If tie, deaults to first max
%     [~,beh] = max([sum(dummy==0),sum(dummy==1),sum(dummy==2)]);
%     % Color upper right square according to behavior
%     frame(1:20,621:640,:) = cat(3,ones(20,20).*col(beh,1),ones(20,20).*col(beh,2),ones(20,20).*col(beh,3));
    % Repeat for prediction
    % Use majority event during frame for color
    dummy = eventPred(startInd:stopInd);
    % If tie, deaults to first max
    [~,beh] = max([sum(dummy==0),sum(dummy==1),sum(dummy==2)]);
    % Color next to above color
    frame(1:20,621:640,:) = cat(3,ones(20,20).*col(beh,1),ones(20,20).*col(beh,2),ones(20,20).*col(beh,3));
    writeVideo(vOut,frame);
    k = k +1;
end
close(vOut)
toc


