function [out] = collateData(sdir,searchStr,vars)
%% Initialization
% Get number of strings going to be used in searchStr
nStr = size(searchStr,1);
% Get number of variables to be concatenated
nVar = numel(vars);
% Checks if searchStr has a criterion per string (stored in second column);
% if not, assumes inclusion and appends column of 'in'
if size(searchStr,2) == 1
    searchStr = [searchStr,repmat({'in'},size(searchStr,1),1)];
end
% Preallocate files
if nStr > 1
    files = cell(nStr,1);
else
    files = cell(1,1);
end
%% Cycle through each searchStr and get file names
for sI = 1:size(searchStr,1)
    [files{sI}] = fileSearch(sdir,searchStr{sI,1},searchStr{sI,2});
end
%% Go through each cell of files and collate data from files within into
% one structure per cell of files
for sI = 1:nStr
    % Get number of files within given cell
    nFile = size(files{sI},2);
    for fI = 1:nFile
        load([sdir,files{sI,fI}]); 
        for vI = 1:nVar
            % If given variable exists, add to that matrix
            if exist(vars{nVar},'var')
                
            end
        end
        clearvars -except sI fI files 
    end
end
