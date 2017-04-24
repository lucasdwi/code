[fileStruct] = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\raw\'},{'.mat'});
files = extractfield(fileStruct{1,1},'name');
cd('C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\raw\')
for fi = 1:size(files,2)
   disp(num2str(fi))
   load(files{1,fi})
   if size(LFPTs.data,1) == 8
       LFPTs.data = LFPTs.data([2,4,6,8],:);
       LFPTs.label = LFPTs.label([2,4,6,8]);
       save(files{1,fi},'ad','adfreq','eventTs','fn','LFPTs','n','pl2','TimeSampEr','ts')
   end
   clearvars -except files fi
end