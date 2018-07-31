[data,~,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\modelData\',{'AC_','in';'AT_','in'},{'pow','coh'},'avg','rel');
%%
x = [data{1};data{2}];
y = [zeros(size(data{1},1),1);ones(size(data{2},1),1)];
allInds = 1:size(x,1);
for ii = 1:size(x,1)/3
    testInds = (3*ii-2):(3*ii);
    trainInds = allInds(~ismember(1:size(x,1),testInds));
    looTrainX{ii} = cat(1,x{trainInds,:});
    looTrainY{ii} = y(trainInds,:);
    looTestX{ii} = cat(1,x{testInds,:});
    looTestY{ii} = y(testInds,:);
end
%%
inds = [];
for ii = 1:4
    for jj = 5:9
        inds = [inds;ii,jj];
    end
end
allInds = 1:size(x,1);
for ii = 1:20
    out = [inds(ii,1)*3-2:inds(ii,1)*3,inds(ii,2)*3-2:inds(ii,2)*3];
    in = allInds(~ismember(allInds,out));
    ltoTrainX{ii} = cat(1,x{in,:});
    ltoTrainY{ii} = y(in,:);
    ltoTestX{ii} = cat(1,x{out,:});
    ltoTestY{ii} = y(out,:);
end
save('C:\Users\Pythia\Documents\GreenLab\data\aaberg\aabergModelDataLto.mat','ltoTrainX','ltoTrainY','ltoTestX','ltoTestY')
%%
for ii = 1:20
   load(['C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\lto\lto',num2str(ii),'.mat']) 
   A(ii) = a;
end
%%
load('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\allDataModel\1.mat')
doubleHist((1-real.err).*100,(1-perm.err).*100,'Main','THC vs. CON','xlab','Accuracy (%)')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\log\log1.mat')
r = cellfun(@(x) x.Rsquared.Adjusted,mdl);
nameVect = names({'rOFC','rIL','rNAc','rPL','lNAc','lPL','lOFC','lIL'},{'d','t','a','b','lg','hg'});
[sR,sInd] = sort(r);
nameSort = nameVect(sInd);
%%
load('C:\Users\Pythia\Documents\GreenLab\data\aaberg\leverModel.mat')
for ii = 3%1:4
   figure
   mdl = fitglm(x(1:3:27,sInd(end-(ii-1))),lever(1:3:27));
   scatter(x(1:3:27,sInd(end-(ii-1))),lever(1:3:27),'ok','filled')
   lsline
   title([nameVect(sInd(end-(ii-1)))])
   xAx = get(gca,'XTick');
   text(xAx(end-1),55,['R^2 = ',num2str(round(mdl.Rsquared.Ordinary,2)),newline,'p = ',num2str(round(table2array(mdl.Coefficients(2,4)),2))])
   xlabel('Normalized Coherence')
   ylabel('Average Lever Presses')
end
%% Grab all power and coherence
files{1} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\modelData\','AC_','in');
files{2} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\aaberg\processed\modelData\','AT_','in');
p = cell(1,2);
c = cell(1,2);
for ii = 1:size(files,2)
   for jj = 1:size(files{ii},2)
       load(files{ii}{jj})
       p{ii} = cat(3,p{ii},psdTrls{1,1}.Overall);
       c{ii} = cat(3,c{ii},mean(coh{1,1}.Cxy,3));
   end
end
% Average
mP = cellfun(@(x) mean(x,3),p,'UniformOutput',0);
mC = cellfun(@(x) mean(x,3),c,'UniformOutput',0);
% STD
sP = cellfun(@(x) std(x,[],3),p,'UniformOutput',0);
sC = cellfun(@(x) std(x,[],3),c,'UniformOutput',0);
% Grab frequency axes
cF = coh{1}.f;
pF = psdTrls{1}.f;
%% Plot rPL-lOFC coherence- with area under curve
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    inds = [nearest_idx3(bands{ii,2}(1),cF) nearest_idx3(bands{ii,2}(2),cF)];
    h = area(cF(inds(1):inds(2)),mC{2}(21,inds(1):inds(2)));
    h.FaceColor = [0.85 0.85 0.85];
    
    h = area(cF(inds(1):inds(2)),mC{1}(21,inds(1):inds(2)));
    h.FaceColor = [0.75 0.75 0.75];
    
    text(sum(bands{ii,2})/2,0.03,latinToGreek({bands{ii,1}}),'HorizontalAlignment','center')
end
clear h
h(1) = plot(cF,mC{1}(21,:),'-r');
h(2) = plot(cF,mC{2}(21,:),'-b');
legend(h(1:2),{'Control','THC'})

title('lNAc-lOFC Coherence')
xlabel('Frequency (Hz)')
ylabel('Magnitude Squared Coherence')
ylim([0 1])
%% Plot lIL power - with area under curve
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    inds = [nearest_idx3(bands{ii,2}(1),pF) nearest_idx3(bands{ii,2}(2),pF)];
    h = area(pF(inds(1):inds(2)),mP{2}(8,inds(1):inds(2)));
    h.FaceColor = [0.85 0.85 0.85];
    
    h = area(pF(inds(1):inds(2)),mP{1}(8,inds(1):inds(2)));
    h.FaceColor = [0.75 0.75 0.75];
    
    text(sum(bands{ii,2})/2,-3,latinToGreek({bands{ii,1}}),'HorizontalAlignment','center')
end
clear h
h(1) = plot(pF,mP{1}(8,:),'-r');
h(2) = plot(pF,mP{2}(8,:),'-b');
legend(h,{'Control','THC'},'Location','sw')
title('lIL Power')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
%% Plot rPL-lOFC coherence - with rectangle
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    rectangle('Position',[bands{ii,2}(1) 0.3 diff(bands{ii,2}) 0.47],'FaceColor',[0.85 0.85 0.85])
end
plot(cF,mC{1}(21,:),'-r')
plot(cF,mC{2}(21,:),'-b')
legend({'Control','THC'})
title('lNAc-lOFC Coherence')
xlabel('Frequency (Hz)')
ylabel('Magnitude Squared Coherence')
ylim([0 1])

%% Plot lIL power - with rectangle
figure
hold on
% Add frequency bands
bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
for ii = 1:size(bands,1)
    rectangle('Position',[bands{ii,2}(1) -85 diff(bands{ii,2}) 50],'FaceColor',[0.85 0.85 0.85])
end
plot(pF,mP{1}(8,:),'-b')
plot(pF,mP{2}(8,:),'-r')
legend({'Control','THC'})
title('lIL Power')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')