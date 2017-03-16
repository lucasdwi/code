function [data,rawPowArr,meanPowArr,relPowArr,rawCohArr,meanCohArr,relCohArr,nTrls,relPowAvg,relCohAvg] = collateData(sdir,searchStr)
% Get file structure for each searchStr
[fileStruct] = fileSearch(sdir,searchStr);
% Load first file to get number of events using size of trls
load([sdir{1},fileStruct{1}(1).name]);
nEvents = size(trls,2);
% Preallocate matrices based on size of fileStruct and nEvents
rawPowArr = cell(1,numel(fileStruct),nEvents);
meanPowArr = cell(numel(fileStruct),nEvents);
relPowArr = cell(numel(fileStruct),nEvents);
rawCohArr = cell(1,numel(fileStruct),nEvents);
meanCohArr = cell(numel(fileStruct),nEvents);
relCohArr = cell(numel(fileStruct),nEvents);
nTrls = cell(numel(fileStruct),nEvents);
%% Go through each file structure and extract data from files therein
for fsi = 1:numel(fileStruct)
    disp(fsi)
    % Preallocate fName
    fName = cell(1,size(fileStruct{1,fsi},1));
    for fi = 1:size(fileStruct{fsi},1)
        disp(fi)
        % Splits filename at '_' and stores first part
        parts = strsplit(fileStruct{fsi}(fi).name,'_');
        fName{fsi,fi} = parts{1};
        % Skip loading first file since it was already loaded for
        % preallocation
        if ~isequal([fsi,fi],[1,1])
            load([sdir{1},fileStruct{fsi}(fi).name]);
        end
        % Convert power structures into cell arrays
        pow{1,1} = struct2cell(psdTrls.event1);
        if isfield(psdTrls,'event2')
            pow{1,2} = struct2cell(psdTrls.event1);
        end
        if isstruct(relPower)
            % Transpose into row vector
            relPower = struct2cell(relPower)';
        end
        % Cycle through all events (usually one or two)
        for ei = 1:size(pow,2)
            % Grab all power spectra (1st cell in Pow) 
            rawPowArr{fsi,fi,ei} = cat(3,pow{1,ei}{1,1}{1,:});
            % Grab mean power spectra
            meanPowArr{fsi,ei} = cat(3,meanPowArr{fsi,ei},pow{1,ei}{2,1});
            % Grab relative power
            relPowArr{fsi,ei} = cat(3,relPowArr{fsi,ei},relPower{1,ei});
            
            % Grab all coherence
            rawCohArr{fsi,fi,ei} = coh{1,ei}.Cxy;  %#ok<USENS>
            % Grab mean coherence
            meanCohArr{fsi,ei} = cat(3,meanCohArr{fsi,ei},coh{1,ei}.mCxy);
            % Grab relative coherence
            relCohArr{fsi,ei} = cat(3,relCohArr{fsi,ei},coh{1,ei}.rel);
            
            % Grab number of trials
            nTrls{fsi,ei} = [nTrls{fsi,ei},size(pow{1,ei}{1,1},2)];
        end
    end
end
%% Replace Infs with NaN
for fsi = 1:numel(fileStruct)
    for ei = 1:nEvents
        meanPowArr{fsi,ei}(isinf(meanPowArr{fsi,ei})) = NaN;
        meanCohArr{fsi,ei}(isinf(meanCohArr{fsi,ei})) = NaN;
    end
end
%% Average relArrs across third dimensions (trials)
relPowAvg = cellfun(@(x) nanmean(x,3),relPowArr,'Un',0);
relCohAvg = cellfun(@(x) nanmean(x,3),relCohArr,'Un',0);
%% Set up data matrices per group: observation (file) X variable
% Preallocate data
data = cell(numel(fileStruct),nEvents);
for fsi = 1:numel(fileStruct)
    for ei = 1:nEvents
        % Reshape and transpose power and coherence matrices to get in
        % right dimensions
        thisPow = reshape(relPowArr{fsi,ei},[size(relPowArr{fsi,ei},1)*size(relPowArr{fsi,ei},2),size(relPowArr{fsi,ei},3)])';
        thisCoh = reshape(relCohArr{fsi,ei},[size(relCohArr{fsi,ei},1)*size(relCohArr{fsi,ei},2),size(relCohArr{fsi,ei},3)])';
        % Concatenate together
        data{fsi,ei} = [thisPow,thisCoh];
    end
end
