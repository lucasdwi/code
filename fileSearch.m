function [fNames] = fileSearch(varargin)
%% Searchs for files in directory given - uses strings given as either inclusion or exclusion criteria
% fileSearch(DIR,STR1,CRIT1,...)
% DIR = directory to search in
% STR1 = string to search for
% CRIT1 = whether to include or exclude files with above string (STR1);
%   either 'in' or 'ex'

% fileSearch('C:\Users\Lucas\Desktop\','foo','in','bar','ex')
% Returns cell array of filenames found on Desktop with 'foo' in the name
% but without 'bar' (i.e. would include file 'foofoo', but not 'foobar')

% fileSearch('C:\Users\Lucas\Desktop\','foo')
% Returns cell arry of all files with 'foo' in name - defaults to inclusion
%% Written by LLD 2016; Influenced by FindFiles by J Lilly (jlab).
%%
% Go to source directory
sdir = varargin{1};
cd(sdir)
% Grab all extensions and criteria - put into vars
vars = varargin(2:end);
%%
% If only one extension and no criterion defined, set to default 'in'
if size(vars,2) == 1
    vars{1,2} = 'in';
end
% Check for proper number of inputs, one criterion per extension
nVar = size(vars,2);
if mod(nVar,2) ~= 0
    error('Incorrect number of inputs! Makes sure each extension has a paired criterion.')
end
% Grab indices of in and ex criteria
inInd = logicFind('in',vars,'==');
exInd = logicFind('ex',vars,'==');
% Use indices to sort vars into proper group
inStr = vars(inInd-1);
exStr = vars(exInd-1);
% Set up fileStruct
inFiles = [];
% Go through strings and populate fileStruct starting with any 'in' strings
if ~isempty(inStr)
    for iSi = 1:size(inStr,2)
        inFiles = [inFiles; dir(['*',inStr{iSi},'*'])]; %#ok<AGROW>
    end
else
    % Otherwise, grab all files
    inFiles = dir('*');
end
% Grab file names
fNamesIn = extractfield(inFiles,'name');
if ~isempty(exStr)
    outFiles = [];
    for eSi = 1:size(exStr,2)
        outFiles = [outFiles; dir(['*',exStr{eSi},'*'])]; %#ok<AGROW>
    end
    fNamesOut = extractfield(outFiles,'name');
end
% If any fNamesOut, remove from fNamesIn
if ~isempty(exStr)
    fNames = fNamesIn(~ismember(fNamesIn,fNamesOut));
else
    fNames = fNamesIn;
end
% Check for, and remove, any non-unique entries
fNames = unique(fNames);
% Check for, and remove, directory entries ('.' or '..')
fNames(logicFind(1,cellfun(@(x) strcmpi('.',x),fNames),'==')) = [];
fNames(logicFind(1,cellfun(@(x) strcmpi('..',x),fNames),'==')) = [];