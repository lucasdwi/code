% fileScratcher: loads all files in scratch
function fileScratcher(sdir)
files = fileSearch(sdir,'.mat');
for ii = 1:numel(files)
    load(files{ii})
end