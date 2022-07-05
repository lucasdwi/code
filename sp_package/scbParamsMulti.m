function [cfg] = scbParamsMulti(file)
%%
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\dianaNVHL\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\splitMat\';
% cfg.sdir = 'D:\paper2\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\maleFemale\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angela\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\aaberg\splitMat\';
% cfg.sdir = 'F:\lsdIRDM\toProcess\';
% cfg.sdir = 'F:\lsdStim\toProcess\';
% cfg.sdir = 'F:\lsd\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\paper3\water\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angela\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\irdm\toProcess\';
% cfg.sdir = 'E:\aaberg\toProcess\';
% cfg.sdir = 'E:\dualSite\toProcess\dualSite\';
% cfg.sdir = 'D:\dualSite\toProcess\singleSite\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\lsdStim\toProcess';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\pulseStim\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\angelaK\toProcess\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\irdm\split\ddt\';
% cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\irdm\toProcess\';
% cfg.sdir = 'G:\GreenLab\data\irdmStim\mat\';
% cfg.sdir = 'G:\GreenLab\data\lsdStim\toProcess\';
% cfg.sdir = 'G:\GreenLab\data\irdmNew\toProcess\';
% cfg.sdir = 'D:\irdm\matNew\';
% cfg.sdir = 'F:\paper3\waterAlcohol\mat\';
% cfg.sdir = 'F:\irdmRound2\toProcess\';
cfg.sdir = 'F:\angelaAlcoholDual\New folder';
% cfg.sdir = 'F:\angelaAlcoholDual\Pre\';
% file = 'H10RegChowNov23';
cfg.file = file;
cfg.nFilt = [57 63];
cfg.dsf = 1;
cfg.thresh = 2; 

cfg.onset = 0.0125;
cfg.offset = 1;
cfg.foi = [1 1 100];
cfg.fixed = 1;
cfg.discrete = 0;
cfg.bands = {'delta',[1,4];
             'theta',[5,10];
             'alpha',[11,14];
             'beta',[15,30];
             'lgamma',[45,65];
             'hgamma',[70,90]};
cfg.overlap = 0.5;
cfg.cohMethod = 'mtm';
cfg.skip = [];
% cfg.eoi = [{'Base1',[0,5]};{'Base2',[0,5]};{'Stim1',[0,5]}];
% for ii = 1:540
%    cfg.eoi(ii,:) =  {['Stim',num2str(ii)],[0 0.25]};
% end
% cfg.eoi(end+1,:) = {'Base',[0 0.25]};
% cfg.eoi = [{'Base',[0,5]};{'Stim',[0,5]}];
% cfg.eoi = [{'Base',[0,5]};{'Washout',[0,5]};{'Stim',[0,5]}];
% cfg.eoi = [{'Base1',[0,5]};{'Base2',[0,5]};{'Base3',[0,5]};...
%     {'Base4',[0,5]};{'Base5',[0,5]};{'Base6',[0,5]};{'Base7',[0,5]};...
%     {'Base8',[0,5]};{'Stim1',[0,5]};{'Stim2',[0,5]};{'Stim3',[0,5]};...
%     {'Stim4',[0,5]};{'Stim5',[0,5]};{'Stim6',[0,5]};{'Stim7',[0,5]}]; 
% cfg.eoi = [{'binge (s'},[-4,1];{'binge (s'},[-3,2];{'binge (s'},[-2,3];...
%     {'binge (s'},[-1,4]];
% cfg.eoi = {'rest',[0 3]};
% cfg.eoi = {'~Both',[0 5]};
% cfg.eoi = {'water',[0 5];'alcohol',[0 5]};
cfg.eoi = {'all',[0 5]};
% cfg.eoi = {'Inter',[0 5]};
% cfg.eoi = {'pre',[0 5];'post',[0 5]};
% cfg.eoi = {'~drink',[0 5]};
% cfg.eoi = {'drink',[0 5];'~drink',[0 5]};
% cfg.eoi = {'trial',[-5 0]};
% for ii = 1:35
%    cfg.eoi(ii,:) = {'Alcohol (S',[-31+ii -26+ii]}; 
% end
% for ii = 1:15
%    cfg.eoi(35+ii,:) = {'Alcohol (E',[-8+ii -3+ii]}; 
% end
% cfg.eoi(end+1,:) = {'~binge',[0 5]};
% cfg.eoi(end+1,:) = {'Binge',[0 5]};
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

% Post 
% for ii = 1:20
%    cfg.eoi(ii,:) = {'binge (e',[ii-10 ii-5]}; 
% end

% 15 seconds pre and post drink and 5 seconds at beginning and end of drink
% for ii = 1:16
%    eoi(ii,:) = {'drink (s',[-ii+1 6-ii]}; 
% end
% for ii = 1:16
%     eoi(ii+16,:) = {'drink (e',[ii-6 ii-1]};
% end
cfg.vis = 'n';
cfg.saveParent = 'F:\angelaAlcoholDual\processedCohort1\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\irdm\processed\ddtMinus1\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\pulseStim\processed\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\paper2\preBingeInfluence\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\dianaNVHL\processed\';
% cfg.saveParent = 'F:\irdmRound2\processedContinuous5\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\angela\processed';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\';
% saveParent = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\processed\';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\paper2\processedAll\';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\dualSite\processed\';
% saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\stimParam\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\processed\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed';
% cfg.saveParent = 'F:\lsdIRDM\processed\lsdImag\';
% cfg.saveParent = 'D:\dualSite\processed\toUseSingleImag\';
% cfg.saveParent = 'F:\lsd\processed\imagWater\';
% cfg.saveParent = 'F:\lsdStim\processed\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\paper3\waterProcessed';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\irdm\processedStim';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\maleFemale\3SecondRest\';
% cfg.saveParent = 'E:\aaberg\processedNew\';
% cfg.saveParent = 'E:\dualSite\processed\';
% cfg.saveParent = 'G:\GreenLab\data\lsdStim\processed\';
% cfg.saveParent = 'C:\Users\Pythia\Documents\GreenLab\data\angelaK\processed\'
% cfg.saveParent = 'G:\GreenLAb\data\irdmNew\processedContinuous\';