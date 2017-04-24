function [] = dbsComp(sdir)
%% Designed to search for 'Pre' and 'Post' files in sdir
% Pre = average and standard deviation power and coherence
preFiles = dir('*Pre*');
for f = 1:size(preFiles,1)
    load(preFiles(f).name)
    prePowAvg(:,:,f) = psdTrls.event1.Overall;
    prePowStd(:,:,f) = psdTrls.event1.OverallStd;
    % Transpose coherence to get in same general format(channel pair x freq)
    preCohAvg(:,:,f) = coh.mCxy';
    preCohStd(:,:,f) = coh.sdCxy';
end
% Post = PSD and coherence for each trial
postFiles = dir('*Post*');
for f = 1:size(postFiles,1)
    load(postFiles(f).name)
    postPow() = ;
    postCoh() = ;
end
%%
for f = 1:length(fType)
    files.(fType{f}) = dir(strcat('*',fType{f},'*'));
end
% Check directory for even number of files
if ~isequal(size(files.(fType{1}),1),size(files.(fType{2}),1))
    error('Warning: There are different number of files between conditions, comparison will be incorrect.')
end
% Load each file and extract power and coherence information
for f = 1:length(fType)
    for fi = 1:size(files.(fType{f}),1)
        load(files.(fType{f})(fi).name)
        disp(['Loading ',files.(fType{f})(fi).name])
    end
end