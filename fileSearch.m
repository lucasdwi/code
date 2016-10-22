function [fileStruct] = fileSearch(sdir,searchStr)
%% Uses search string(s) to grab filenames from directory
% INPUTS
% sdir = source directory; format = cell
% searchStr = string to use to select files; format = cell
%%
% Go to directory
fileStruct = cell(1,numel(searchStr));
for sdi = 1:numel(sdir)
    cd(sdir{sdi})
    for ssi = 1:numel(searchStr)
        fileStruct{ssi} = dir(['*',searchStr{ssi},'*']);
    end
end
