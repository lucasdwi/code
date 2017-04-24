function [fileStruct] = fileSearch(sdir,searchStr)
%% Uses search string(s) to grab filenames from directory
% INPUTS
% sdir = source directory; format = cell
% searchStr = string to use to select files; format = cell
%%
fileStruct = cell(1,numel(searchStr));
for sdi = 1:numel(sdir)
    % Go to directory
    cd(sdir{sdi})
    for ssi = 1:numel(searchStr)
        % Get file details
        fileStruct{ssi} = dir(['*',searchStr{ssi},'*']);
    end
end
% % If only one cell, get rid of singleton dimension
% if numel(searchStr) == 1
%     fileStruct = fileStruct{1};
% end
