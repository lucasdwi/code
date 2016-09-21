function [] = dbsComp(sdir)
%% Cd to file directory
cd(sdir)
%% Designed to search for 'Pre' and 'Post' files in sdir
% Pre = average and standard deviation power and coherence
preFiles = dir('*Pre*');
for f = 1:size(preFiles,1)
    load(preFiles(f).name)
    % Dim = channel X freq X file
    prePowAvg(:,:,f) = psdTrls.event1.Overall;
    prePowStd(:,:,f) = psdTrls.event1.OverallStd;
    % Dim = channel pair X freq X file
    preCohAvg(:,:,f) = coh.mCxy;
    preCohStd(:,:,f) = coh.sdCxy;
    % Extract frequency axis from last file for plotting
    if f == size(preFiles,1)
        powFreq1 = psdTrls.F;
    end
end
% Post = PSD and coherence for each trial
postFiles = dir('*Post*');
for f = 1:size(postFiles,1)
    load(postFiles(f).name)
    % Dim = channel X freq X trial
    postPow{f}(:,:,:) = cat(3,psdTrls.event1.Pow{1,:});
    % Dim = channel pair X freq X trial
    postCoh{f}(:,:,:) = coh.Cxy;
    % Extract frequency axis from last file for plotting
    if f == size(postFiles,1)
        powFreq2 = psdTrls.F;
    end
end
% Check for same number of files in both pre and post conditions
if ~isequal(size(preFiles,1),size(postFiles,1))
    error('There are not the same number of files in the two conditions, comparison will not work.')
end
if ~isequal(powFre1,powFreq2)
    error('There are different frequency axes for power in the two conditions, files not comparable.')
end
%% Pre to post comparison normalized by SD
for f = 1:size(postFiles,1)
   for t = 1:size(postPow{f},3)
       % Dim = channel(pair) X freq X trial
       powChange{f}(:,:,t) = (sq(postPow{f}(:,:,t)) - sq(prePowAvg(:,:,f)))./prePowStd(:,:,f);
       cohChange{f}(:,:,t) = (sq(postCoh{f}(:,:,t)) - sq(preCohAvg(:,:,f)))./preCohStd(:,:,f);
   end
end
%% Plot change in power and coherence


for f = 1:size(postFiles,1)
    figure; hold on;
    for c = 1:size(powChange{1},1)
        subplot(2,2,c)
        imagesc(1:size(powChange{f},3),powFreq1,sq(powChange{f}(c,:,:)));
    end
end
%%
% for f = 1:length(fType)
%     files.(fType{f}) = dir(strcat('*',fType{f},'*'));
% end
% % Check directory for even number of files
% if ~isequal(size(files.(fType{1}),1),size(files.(fType{2}),1))
%     error('Warning: There are different number of files between conditions, comparison will be incorrect.')
% end
% % Load each file and extract power and coherence information
% for f = 1:length(fType)
%     for fi = 1:size(files.(fType{f}),1)
%         load(files.(fType{f})(fi).name)
%         disp(['Loading ',files.(fType{f})(fi).name])
%     end
% end