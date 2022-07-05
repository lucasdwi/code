%% MIA with alcohol injection
% convert to matlab files
ConvertPl2All_Files('C:\Users\GreenLab\Desktop\angela\MIA\MIAwAlcInj\PRE');
% preprocess pre inj data
fileCycle('scb',{'.mat'},[],'C:\Users\GreenLab\Desktop\angela\MIA\MIAwAlcInj\PRE\mat');
% preprocess post inj data
fileCycle('scb',{'.mat'},[],'C:\Users\GreenLab\Desktop\angela\MIA\MIAwAlcInj\POST\mat');