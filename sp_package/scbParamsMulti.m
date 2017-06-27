function [sdir,file,filter,dsf,thresh,onset,offset,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent] = scbParamsMulti(file)
%%
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\paper2\toProcess\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\angela\mat\nicks\';
%sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\baseData\';
%sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\';
% sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\raw\';
%file = 'H10RegChowNov23';
file = file;
filter = 'y';
dsf = 5;
thresh = 2; 
onset = 0.0125;
offset = 40;
foi = [1 2 100];
% bands = {'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
cycles = 3;
ftimwin = [];
overlap = 0.5;
cohMethod = 'mat';
eoi = {'binge',[0 5];'rest',[0 5]};
% for ii = 1:61
%    eoi(ii,:) = {'binge (s',[-4-ii 1-ii]}; 
% end
% eoi = {'binge (s',[-60 -55]};%{'binge (s',[-11 -6];'binge (s',[-12 -7];'binge (s',[-13 -8];'binge (s',[-14 -9];'binge (s',[-15 -10];'binge (s',[-16 -11];'binge (s',[-17 -12];'binge (s',[-18 -13];'binge (s',[-19 -14];'binge (s',[-20 -15];'binge (s',[-21 -16];'binge (s',[-22 -17]};
saveParent = 'C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\';