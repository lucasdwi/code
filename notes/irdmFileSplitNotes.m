oldLFPTs = LFPTs; oldEventTs = eventTs;
%%
LFPTs.data = oldLFPTs.data(1:8,:); 
LFPTs.label = oldLFPTs.label(1:8);
eventTs.t = oldEventTs.t(25:32);
eventTs.label = oldEventTs.label(25:32);
% eventTs.t = oldEventTs.t(33:40);
% eventTs.label = oldEventTs.label(33:40);
save('IRDM18_DDT_2019-04-09.mat','LFPTs','eventTs','pl2','adfreq')
%%
LFPTs.data = oldLFPTs.data(9:16,:);
LFPTs.label = oldLFPTs.label(9:16);
% eventTs.t = oldEventTs.t(49:56);
% eventTs.label = oldEventTs.label(49:56);
eventTs.t = oldEventTs.t(33:40);
eventTs.label = oldEventTs.label(33:40);
save('IRDM11_DDT_2019-04-09.mat','LFPTs','eventTs','pl2','adfreq')
%%
LFPTs.data = oldLFPTs.data(17:24,:);
LFPTs.label = oldLFPTs.label(17:24);
% eventTs.t = oldEventTs.t(41:48);
% eventTs.label = oldEventTs.label(41:48);
eventTs.t = oldEventTs.t(49:56);
eventTs.label = oldEventTs.label(49:56);
save('IRDM24_DDT_2019-04-09.mat','LFPTs','eventTs','pl2','adfreq')
%%
LFPTs.data = oldLFPTs.data(25:32,:);
LFPTs.label = oldLFPTs.label(25:32);
eventTs.t = oldEventTs.t(49:56);
eventTs.label = oldEventTs.label(49:56);
save('IRDM21_DDT_2019-04-09.mat','LFPTs','eventTs','pl2','adfreq')
%%
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\irdm\split\ddt\','.mat');
for ii = 1:size(files,2)
    load(files{ii})
    adfreq = 2000;
    save(files{ii},'LFPTs','eventTs','pl2','adfreq')
end