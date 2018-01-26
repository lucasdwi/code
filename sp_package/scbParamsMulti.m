function [sdir,file,filter,dsf,thresh,onset,offset,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent] = scbParamsMulti(file)
%%
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\mike\toProcess\';
sdir = 'C:\Users\Pythia\Documents\GreenLab\data\paper2\toProcess\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angela\mat\nicks\';
%sdir = 'C:\Users\Pythia\Documents\GreenLab\data\megan\baseData\';
%sdir = 'C:\Users\Pythia\Documents\GreenLab\data\WilderBinge\channel_renamed\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\twoSiteStim\mat\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\maleFemale\mat\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angela\toProcess\';
%file = 'H10RegChowNov23';
file = file;
filter = 'y';
dsf = 5;
thresh = 2; 
onset = 0.0125;
offset = 40;
foi = [1 2 200];
% bands = {'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
% bands = {'lgam',[45,65]};
cycles = 3;
ftimwin = [];
overlap = 0.5;
cohMethod = 'mat';
eoi = {'all',[0 5]};
% eoi = {'rest',[0 5]};
% eoi = {'Base',[0 5];'Int1',[0 5];'Int2',[0 5];'Int3',[0 5];'Int4',[0 5];'Int5',[0 5];'Post',[0 5]};
% eoi = {'binge',[0 5];'notbinge',[0 5]};
% for ii = 1:61
%    eoi(ii,:) = {'binge (s',[-4-ii 1-ii]}; 
% end
% for ii = 1:31
%    eoi(ii,:) = {'binge (s',[-1+ii 4+ii]}; 
% end
% for ii = 1:61
%    eoi(ii,:) = {'binge (e',[-10+ii -5+ii]}; 
% end
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\angela\processed';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\';
% saveParent = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\processed\';
saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\paper2\processedAll\';