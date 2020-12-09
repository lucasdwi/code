%% Compare loading times between local ssd, local hd, and network drive
%% Local ssd
tic
files = filSearch('H:\Shared drives\Green Lab2\testFiles\','.mat');
for ii = 1:numel(files)
    load(files{ii})
    ssd(ii) = toc;
end
%% Local hd
files = filSearch('F:\testFiles\','.mat');
%% Network drive
files = filSearch('E:\testFiles\','.mat');