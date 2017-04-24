function [powChange,cohChange] = dbsComp(sdir)
%% Cd to file directory
cd(sdir)
%% Designed to search for 'Pre' and 'Post' files in sdir
% Pre = average and standard deviation power and coherence
tic
preFiles = dir('*Pre*');
for f = 1:size(preFiles,1)
    load(preFiles(f).name)
    % Use first file to preallocate and extract freq axis; as long as all
    % files are from the same processing pipeline all structures/size
    % should match
    if f == 1
        prePowAvg = zeros(horzcat(size(psdTrls.event1.Overall),size(preFiles,1)));
        prePowStd = zeros(size(prePowAvg));
        preCohAvg = zeros(horzcat(size(coh.mCxy),size(preFiles,1)));
        preCohStd = zeros(size(preCohAvg));
        powFreq1 = psdTrls.F;
        cohFreq1 = coh.freq;
    end
    % Dim = channel X freq X file
    prePowAvg(:,:,f) = psdTrls.event1.Overall;
    prePowStd(:,:,f) = psdTrls.event1.OverallStd;
    % Dim = channel pair X freq X file
    preCohAvg(:,:,f) = coh.mCxy;
    preCohStd(:,:,f) = coh.sdCxy;
end
% Post = PSD and coherence for each trial
postFiles = dir('*Post*');
for f = 1:size(postFiles,1)
    load(postFiles(f).name)
    % Use first file to preallocate and extract freq axis; as long as all
    % files are from the same processing pipeline all structures/size
    % should match
    if f == 1
        postPow = cell(size(postFiles,1));
        postCoh = cell(size(postFiles,1));
        postTrialTimes = cell(1,size(postFiles,1));
        chk_nan = cell(size(postFiles,1));
        powFreq2 = psdTrls.F;
        cohFreq2 = coh.freq;
    end
    % Dim = channel X freq X trial
    postPow{f}(:,:,:) = cat(3,psdTrls.event1.Pow{1,:});
    % Dim = channel pair X freq X trial
    postCoh{f}(:,:,:) = coh.Cxy;
    % Extract sampleInfo for time axis
    postTrialTimes{f}(:,:) = trls{1}.sampleinfo;
    % Extract chk_nan for plotting blank space
    chk_nan{f} = hist.chk_nan;
    % Extract labels
    labels{f} = trls{1}.label;
    % Extract frequency axis from last file for plotting
    if f == size(postFiles,1)
        powFreq2 = psdTrls.F;
        fs = hist.adfreq;
        minInt = hist.minInt;
    end
end
% Check for same number of files in both pre and post conditions
if ~isequal(size(preFiles,1),size(postFiles,1))
    error('There are not the same number of files in the two conditions, comparison will not work.')
end
if ~isequal(powFreq1,powFreq2)
    error('Power frequency axes are different between conditions.')
end
if ~isequal(cohFreq1,cohFreq2)
    error('Coherence frequency axes are different between conditions.')
end
toc
%% Pre to post comparison normalized by SD
% Preallocate
powChange = cell(1,size(postFiles,1));
cohChange = cell(1,size(postFiles,1));
for f = 1:size(postFiles,1)
   for t = 1:size(postPow{f},3)
       % Dim = channel(pair) X freq X trial
       powChange{f}(:,:,t) = (postPow{f}(:,:,t) - prePowAvg(:,:,f))./prePowStd(:,:,f);
       cohChange{f}(:,:,t) = (postCoh{f}(:,:,t) - preCohAvg(:,:,f))./preCohStd(:,:,f);
   end
end
%% Check for incomplete postTrialTimes
for f = 1:size(postFiles,1)
    dummyP{f} = []; dummyC{f} = [];
    %if chk_nan{f} ~= 0
    % Find where gaps are
    diffInd = find(diff(postTrialTimes{f}(:,1))~= minInt*fs);
    % Copy up to first gap
    if ~isempty(diffInd)
        newTime{f} = postTrialTimes{f}(1:diffInd(1),:);
        for di = 1:length(diffInd)
            newTime{f}(end+1,:) = [newTime{f}(end,2)+1,postTrialTimes{f}(diffInd(di)+1,1)-1];
            % Copy up to next gap unless at last gap already, then copy to
            % end
            if di ~= length(diffInd)
                newTime{f} = vertcat(newTime{f},postTrialTimes{f}(diffInd(di)+1:diffInd(di+1),:));
            else if di == length(diffInd)
                    newTime{f} = vertcat(newTime{f},postTrialTimes{f}(diffInd(di)+1:end,:));
                end
            end
        end
        % Check if newTime starts at 1
        if newTime{f}(1,1) ~= 1
            newTime{f} = vertcat([1,newTime{f}(1,1)-1],newTime{f});
            % Add index to diffInd last so that it is last to be change in
            % data
            addOneInd = 'y';
        else addOneInd = 'n';
        end
        % Check that newTime intervals increase by one each step
        for ii = 1:size(newTime{f},1)-1
            testStep(ii) = newTime{f}(ii+1,1)-newTime{f}(ii,2);
        end
        if sum(testStep)/length(testStep) ~= 1
            error('Intervals are not increasing correctly.')
        end
        % Add NaN pages to powChange and cohChange corresponding to diffInd
        % Uses dummy variables as placeholders until complete, then
        % overwrites old powChange and cohChange
        for di = 1:length(diffInd)
            if di == 1
                dummyP{f} = cat(3,powChange{f}(:,:,1:diffInd(di)),NaN(size(powChange{f},1),size(powChange{f},2)));
                dummyC{f} = cat(3,cohChange{f}(:,:,1:diffInd(di)),NaN(size(cohChange{f},1),size(cohChange{f},2)));
            else
                if di == length(diffInd)
                    % If last, go to end
                    dummyP{f} = cat(3,dummyP{f},powChange{f}(:,:,diffInd(di-1)+1:diffInd(di)),NaN(size(powChange{f},1),size(powChange{f},2)),powChange{f}(:,:,diffInd(di)+1:end));
                    dummyC{f} = cat(3,dummyC{f},cohChange{f}(:,:,diffInd(di-1)+1:diffInd(di)),NaN(size(cohChange{f},1),size(cohChange{f},2)),cohChange{f}(:,:,diffInd(di)+1:end));
                else
                    dummyP{f} = cat(3,dummyP{f},powChange{f}(:,:,diffInd(di-1)+1:diffInd(di)),NaN(size(powChange{f},1),size(powChange{f},2)));
                    dummyC{f} = cat(3,dummyC{f},cohChange{f}(:,:,diffInd(di-1)+1:diffInd(di)),NaN(size(cohChange{f},1),size(cohChange{f},2)));
                end
            end
        end
        % If added one index, then add that NaN page
        if addOneInd == 'y'
            dummyP{f} = cat(3,NaN(size(powChange{f},1),size(powChange{f},2)),dummyP{f});
            dummyC{f} = cat(3,NaN(size(cohChange{f},1),size(cohChange{f},2)),dummyC{f});
        end
    end
    %end
    % Replace original powChange and cohChange arrays with dummyP and
    % dummyC
    if ~isempty(dummyP{f}) && ~isempty(dummyC{f})
        powChange{f} = dummyP{f};
        cohChange{f} = dummyC{f};
        postTrialTimes{f} = newTime{f};
    end
end
%% Plot change in power
for f = 1:2%size(postFiles,1)
    figure; hold on;
    for c = 1:size(powChange{1},1)
        subplot(2,2,c)
        uimagesc(postTrialTimes{f}(:,1)./(60*fs),powFreq1,sq(powChange{f}(c,:,:)));
        ylim([1 100])
        xlabel('Time (min)'); ylabel('Frequency (Hz)');
        title(['Change in ',labels{f}{c},' Power'])
    end
    % Get name
    thisName = postFiles(f,1).name;
    [newName,remain] = strtok(thisName,'_'); newName = [newName,strtok(remain,'_')];
    mtit(newName)
end
%% Plot change in coherence
% Get channel combinations
cmb = nchoosek(1:4,2);
for f = 1:2%size(postFiles,1)
    for c = 1:size(cmb,1)
    channelCmb(c,:) = labels{f}(cmb(c,:));
    end
    figure; hold on;
    for c = 1:size(cohChange{1},1)
        subplot(3,2,c)
        uimagesc(postTrialTimes{f}(:,1)./(60*fs),cohFreq1,sq(cohChange{f}(c,:,:)));
        ylim([1 100])
        xlabel('Time (min)'); ylabel('Frequency (Hz)');
        title([channelCmb{c,1},' - ',channelCmb{c,2},' Coherence'])
    end
    % Get name
    thisName = postFiles(f,1).name;
    [newName,remain] = strtok(thisName,'_'); newName = [newName,strtok(remain,'_')];
    mtit(newName)
end
