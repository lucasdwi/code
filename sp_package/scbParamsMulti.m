function [sdir,file,nFilt,dsf,thresh,onset,offset,foi,bands,overlap,cohMethod,eoi,vis,saveParent] = scbParamsMulti(file)
%%
sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\splitMat\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\stimParam\toProcess\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\dualSite\toProcess\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\mike\toProcess\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\paper2\toProcess\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angela\mat\nicks\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\megan\baseData\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\WilderBinge\channel_renamed\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\twoSiteStim\mat\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\maleFemale\mat\';
% sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angela\toProcess\';
% file = 'H10RegChowNov23';
file = file;
nFilt = [57 63];
dsf = 5;
thresh = 1; 
onset = 1;
offset = 1;
foi = [1 1 100];
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
overlap = 0.5;
cohMethod = 'mtm';
eoi = {'all',[0 5]};
% eoi = [];
% for ii = 1:2
%     eoi = [eoi;{['Base',num2str(ii),' '],[0 5]}];
% end
% for ii = 1
%    eoi = [eoi;{['Inter',num2str(ii),' '],[0,5]}]; 
% end
% eoi = {'base',[0 5];'post',[0 5]};
% eoi = {'drink',[0 5];'~drink',[0 5]};
% eoi = {'rest',[0 5]};
% eoi = {'Base',[0 5];'Int1',[0 5];'Int2',[0 5];'Int3',[0 5];'Int4',[0 5];'Int5',[0 5];'Post',[0 5]};
% eoi = {'binge',[0 5];'notbinge',[0 5]};

% Pre Drinking
% for ii = 1:61
%    eoi(ii,:) = {'drink (s',[-4-ii 1-ii]}; 
% end

% Post Drinking
% for ii = 1:61
%    eoi(ii,:) = {'drink (e',[ii-1 ii+4]}; 
% end

% 15 seconds pre and post drink and 5 seconds at beginning and end of drink
% for ii = 1:16
%    eoi(ii,:) = {'drink (s',[-ii+1 6-ii]}; 
% end
% for ii = 1:16
%     eoi(ii+16,:) = {'drink (e',[ii-6 ii-1]};
% end
vis = 'n';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\angela\processed';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\';
% saveParent = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\processed\';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\paper2\processedAll\';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\dualSite\processed\';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\stimParam\';
saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\processed\';