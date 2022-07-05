load('F:\irdmRound2\ddtFiles2.mat')
% Remove files with: 'skip','neural','seizure'
skipInd = contains(processedData(:,6),'skip') | ...
    contains(processedData(:,6),'neural') | ...
    contains(processedData(:,6),'seizure') | ...
    contains(processedData(:,6),'OFC');
processedData(skipInd,:) = [];

% Get list of all processed files
files = fileSearch('F:/irdmRound2/processedContinuous3/','.mat');
% Set n, numer of baselines to use
n = 3;
% Get unique IDs
uID = unique(processedData(:,1));
% Go through each animal and find the recordings with interventions and
% their prior baselines; also calculate number of clean trials in data
inter = {'sIL','sNAcC','guan0.3','mph3','mph0.3','mph1','mph0.1',...
    'sNAcC-LSD','sIL-LSD'};
% Clear table
interTable = array2table(zeros(0,10));
% Preallocate
interBaseInd = cell(numel(uID),numel(inter));
lastBase = cell(numel(uID),n);
baseInd = cell(numel(uID),1);
intInd = cell(numel(uID),numel(inter));
for ii = 1:numel(uID)
    disp(ii)
    % Get logical for animal ID matches
    idLog = contains(processedData(:,1),uID{ii});
    idInd = logicFind(1,idLog,'==');
    % Sort by date within each animal
    [theseDates,sortInd] = sort(datetime(processedData(idLog,2)));
    processedData(idLog,:) = processedData(idInd(sortInd),:);
    % Base inds
    baseInd{ii} = logicFind(1,contains(processedData(:,6),...
        'base') & idLog,'==');
    % Find recordings from each animal (ii) with each intervention (jj)
    for jj = 1:numel(inter)
        % For first two interventions, omit entries with LSD
        if jj <= 2
            intInd{ii,jj} = logicFind(1,contains(processedData(:,6),...
                inter{jj}) & idLog & ...
                ~contains(processedData(:,6),'LSD'),'==');
        else
            intInd{ii,jj} = logicFind(1,contains(processedData(:,6),...
                inter{jj}) & idLog,'==');
        end
    end
    % Figure out which is the first intervention by using indices (sorted
    % by date)
    firstInt = min(cell2mat(cellfun(@min,intInd(ii,:),...
        'UniformOutput',false)));
    if ~isempty(firstInt)
        % Figure out last three baselines before any intervention
        preBase = baseInd{ii}(baseInd{ii}<firstInt);
        lastBaseInd = preBase(1,end-(n-1):end);
        % Load last baselines (normalized data)
        for k = 1:n
            thisFile = files(contains(files,...
                processedData{lastBaseInd(k),1}) & ...
                contains(files,processedData{lastBaseInd(k),2}));
            delay{ii,k} = processedData{lastBaseInd(k),9};
            percentDelay{ii,k} = processedData{lastBaseInd(k),10};
            % Load and remove NaNed data
            load(['F:/irdmRound2/processedContinuous3/',thisFile{1}]);
            thisNotNan = ~isnan(squeeze(psdTrls{1}.relPow(1,1,1,:))) &...
                ~isnan(squeeze(coh{1}.normBandCoh(1,1,1,:)));
            notNanPow = squeeze(psdTrls{1}.relPow(:,:,1,thisNotNan));
            notNanCoh = squeeze(coh{1}.normBandCoh(:,:,1,thisNotNan));
            % Fill in NaNs from missing channels if any
            if size(LFPTs.data,1) < 8
                [newNotNanPow,newNotNanCoh] = addNaNs(LFPTs,...
                    notNanPow,notNanCoh,8);
            else
                newNotNanPow = notNanPow;
                newNotNanCoh = notNanCoh;
            end
            % Reshape
            rePow = reshape(newNotNanPow,48,sum(thisNotNan))';
            reCoh = reshape(permute(newNotNanCoh,[2,1,3]),168,...
                sum(thisNotNan))';
            lastBase{ii,k} = [rePow,reCoh];
        end
    end
end
%%
save('F:/irdmRound2/lastBaseNorm2.mat','lastBase','uID','percentDelay','delay')