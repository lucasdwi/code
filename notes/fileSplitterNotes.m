% Step 1: get all pl2 files into a folder
% Step 2: run ConvertPl2All_Files on that folder
ConvertPl2All_Files('F:\DDT\pl2\');
% Step 3: run fileSplitter in that folder; defining your source directory,
% search terms (filter to grab files), and number of channels for each
% animal
fileSplitter('F:\DDT\pl2\','.mat',8)

