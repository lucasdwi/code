% Grab structure of all files
files = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper2\lasso\'},{'.mat'});
% Get list all of animals and find unique
for fi = 1:size(files,1)
    parts = regsplit(files(fi).name);
    animals{fi} = parts{1};
end
uniqueAnimals = unique(animals);
% Use unique animals to go through each animal's set of files
for ui = 1:numel(uniqueAnimals)
   animalFiles = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper2\lasso\'},uniqueAnimals(ui));
   % Find first file
   for ii = 1:size(animalFiles,1)
       name = strsplit(animalFiles(ii).name,'_');
       parts = regsplit(name{1});
       date(ii) = datetime(datestr(parts{end}));
   end
   [~,firstI] = min(date);   
   % Run logistic modeling on all files
   for ii = 1:size(animalFiles,1)
       % Get T1 index of this file
       
   % If first file, keep beta values and indices for future models
   
   % Use first file model to predict other files
end

