<<<<<<< Updated upstream
function [sdir,file,filter,dsf,thresh,onset,offset,minInt,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent] = scbParamsMulti(file)
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\baseData\';
=======
function [sdir,file,filter,dsf,thresh,onset,offset,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent] = scbParamsMulti(file)
%%
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\paper2\mat\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\angela\mat\';
%sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\baseData\';
>>>>>>> Stashed changes
%sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\raw\';
%file = 'H10RegChowNov23';
file = file;
filter = 'y';
dsf = 5;
thresh = 2; 
<<<<<<< Updated upstream
onset = 5;
offset = 17000;
minInt = 5;
foi = [1 2 100];
bands = {'theta',[4,8];'alpha',[8,12];'beta',[13,30];'lgam',[45,65];'hgam',[70,90]};
=======
onset = 0.0125;
offset = 40;
% minInt = 5;
foi = [1 2 100];
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
>>>>>>> Stashed changes
cycles = 3;
ftimwin = [];
overlap = 0.5;
cohMethod = 'mat';
<<<<<<< Updated upstream
eoi = {'rest',[0 5]};
saveParent = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\secondPass\';
=======
eoi = {'binge',[0 5];'notbinge',[0 5]};
saveParent = 'C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed4\';
>>>>>>> Stashed changes
