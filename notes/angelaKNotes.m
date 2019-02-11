files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\angelaK\mat\','.mat');
%% Manually go through each file (fI) 
fI = 10;
load(files{fI})
% Grab this filename
thisFile = files{fI}
chan = 8;
%% Set number of animals based on file name and split
n = 2;
% Remove file extension by finding period; '.'
thisFile = thisFile(1:strfind(thisFile,'.')-1);
% Split str at all underscores; '_'
parts = strsplit(thisFile,'_');

% Store LFPTs
oldLFPTs = LFPTs;
% Go through each sub-file, pull out LFPTs data and save
for ii = 1:n
    clear LFPTs;
    LFPTs.type = oldLFPTs.type;
    LFPTs.tvec = oldLFPTs.tvec;
    LFPTs.data = oldLFPTs.data(chan*ii-chan+1:chan*ii,:);
    LFPTs.label = oldLFPTs.label(chan*ii-chan+1:chan*ii);
    LFPTs.cfg = oldLFPTs.cfg;
    save(strjoin([parts(1),parts(1+ii),parts(end-1:end)],'_'),'LFPTs',...
        'adfreq','eventTs','pl2')
end
%% Collate and run t-tests
[dep5,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\angelaK\processed\'],{'_5_Dep'},{'pow','coh'},'trl','rel');
[nonDep5,samp,files] = collateData(['C:\Users\Pythia\Documents\'...
    'GreenLab\data\angelaK\processed\'],{'_5_NonDep'},{'pow','coh'},...
    'trl','rel');

[dep9,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\angelaK\processed\'],{'_NO-CHANNEL'},{'pow','coh'},'trl','rel');
[nonDep9,samp,files] = collateData(['C:\Users\Pythia\Documents\'...
    'GreenLab\data\angelaK\processed\'],{'_9_NonDep'},{'pow','coh'},...
    'trl','rel');
% Get variable names - account for lost first channel in dep9 recordings
nameVect = names({'rOFC','rIL','rNAc','rPL','lNAc','lPL','lOFC','lIL'},...
    {'d','t','a','b','lg','hg'});
nameVect2 = names({'rIL','rNAc','rPL','lNAc','lPL','lOFC','lIL'},...
    {'d','t','a','b','lg','hg'});
% Get indices of nameVect that correspond to nameVect2
subNameInd = ismember(nameVect,nameVect2);
% Collapse data into single matrices
dep5All = cat(1,dep5{1}{:});
nonDep5All = cat(1,nonDep5{1}{:});
dep9All = cat(1,dep9{1}{:});
nonDep9All = cat(1,nonDep9{1}{:});
% Run t-tests for animal 5
[~,p5,pAdj5] = bulkT(dep5All,nonDep5All,0,'fdr');
[~,p9,pAdj9] = bulkT(dep9All,nonDep9All(:,subNameInd(1:216)),0,'fdr');
% Find significance mismatches
misInd = logicFind(1,(pAdj5(subNameInd(1:216))<=0.05) ~= (pAdj9<=0.05),'==');
%% Make tables
table5 = table(mean(dep5All,1)',std(dep5All,[],1)',mean(nonDep5All,1)',...
    std(nonDep5All,[],1)',pAdj5','RowNames',nameVect(1:216)',...
    'VariableNames',{'DepMean','DepStd','NonDepMean','NonDepStd','pAdj'});

table9 = table(mean(dep9All,1)',std(dep9All,[],1)',...
    mean(nonDep9All(:,subNameInd(1:216)),1)',...
    std(nonDep9All(:,subNameInd(1:216)),[],1)',pAdj9',...
    'RowNames',nameVect2(1:168)',...
    'VariableNames',{'DepMean','DepStd','NonDepMean','NonDepStd','pAdj'});
%%
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\angelaK\processed','5_Dep');
for ii = 1:size(files,2)
    load(files{ii})
    allDepCoh5{ii} = coh{1}.Cxy;
end
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\angelaK\processed','5_Non');
for ii = 1:size(files,2)
    load(files{ii})
    allNonDepCoh5{ii} = coh{1}.Cxy;
end

depCoh5=[];
for ii = 1:3
    depCoh5 = cat(3,depCoh5,allDepCoh5{ii}(14,:,:));
end
nonDepCoh5=[];
for ii = 1:4
    nonDepCoh5 = cat(3,nonDepCoh5,allNonDepCoh5{ii}(14,:,:));
end

figure
hold on
sd = plot(coh{1}.f,mean(depCoh5,3),'b');
sn = plot(coh{1}.f,mean(nonDepCoh5,3),'r');
% sd = shadedErrorBar(coh{1}.f,mean(depCoh,3),std(depCoh,[],3),'-b',1);
% sn = shadedErrorBar(coh{1}.f,mean(nonDepCoh,3),std(nonDepCoh,[],3),'-r',1);
plot([70 70],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([90 90],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([15 15],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([30 30],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
xlabel('Frequency (Hz)'); ylabel('Coherence'); title('AH5')
% legend([sd.mainLine sn.mainLine],'Dep','Non')
legend([sd sn],'Dep','Non')
%%
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\angelaK\processed','9_Dep');
for ii = 1:size(files,2)
    load(files{ii})
    allDepCoh9{ii} = coh{1}.Cxy;
end
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\angelaK\processed','9_Non');
for ii = 1:size(files,2)
    load(files{ii})
    allNonDepCoh9{ii} = coh{1}.Cxy;
end

depCoh9=[];
for ii = 1:3
    depCoh9 = cat(3,depCoh9,allDepCoh9{ii}(7,:,:));
end
nonDepCoh9=[];
for ii = 1:4
    nonDepCoh9 = cat(3,nonDepCoh9,allNonDepCoh9{ii}(7,:,:));
end

figure
hold on
sd = plot(coh{1}.f,mean(depCoh9,3),'b');
sn = plot(coh{1}.f,mean(nonDepCoh9,3),'r');
% sd = shadedErrorBar(coh{1}.f,mean(depCoh,3),std(depCoh,[],3),'-b',1);
% sn = shadedErrorBar(coh{1}.f,mean(nonDepCoh,3),std(nonDepCoh,[],3),'-r',1);
plot([70 70],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([90 90],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([15 15],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([30 30],[0.4 0.9],'--','col',[0.5 0.5 0.5],'LineWidth',1)
xlabel('Frequency (Hz)'); ylabel('Coherence'); title('AH9')
% legend([sd.mainLine sn.mainLine],'Dep','Non')
legend([sd sn],'Dep','Non')
%%
% Plot difference from baseline (nonDepCoh)
figure
h1 = plot(coh{1}.f,mean(depCoh5,3)-mean(nonDepCoh5,3));
hold on
h2 = plot(coh{1}.f,mean(depCoh9,3)-mean(nonDepCoh9,3));
plot([70 70],[-0.15 0.2],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([90 90],[-0.15 0.2],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([15 15],[-0.15 0.2],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([30 30],[-0.15 0.2],'--','col',[0.5 0.5 0.5],'LineWidth',1)
box off
legend([h1 h2],{'Rat1','Rat2'})
title('dep-nonDep')
% Plot percent change from baseline
figure
h1 = plot(coh{1}.f,(mean(depCoh5,3)-mean(nonDepCoh5,3))./...
    mean(nonDepCoh5,3).*100);
hold on
h2 = plot(coh{1}.f,(mean(depCoh9,3)-mean(nonDepCoh9,3))./...
    mean(nonDepCoh9,3).*100);
plot([70 70],[-20 30],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([90 90],[-20 30],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([15 15],[-20 30],'--','col',[0.5 0.5 0.5],'LineWidth',1)
plot([30 30],[-20 30],'--','col',[0.5 0.5 0.5],'LineWidth',1)
box off
legend([h1 h2],{'Rat1','Rat2'})
%title('(dep-nonDep)/nonDep')
title('Percent change from non-dependence')
xlabel('Frequency (Hz)')
ylabel('% Change')