function preProcDBS(fDir)
cd(fDir);
fs = dir('*.mat');
files = extractfield(fs,'name')';
%%
for f = 1:length(files)
    tic
    disp(strcat('Processing file ',files(f)))
    load(files{f})
    LFPTs.data = LFPTs.data(2:2:8,:);
    LFPTs.label = LFPTs.label(2:2:8);
    save(strcat(fDir,'edited\',files{f}));
    toc
end
