%% Megan Notes
%% 36
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\36\';
searchStr = {'base1','in';'base2','in';'alc1','in';'cloz1','in'};
% searchStr = {'base2','in';'alc2','in';'cloz4','in'};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1},searchStr{iS,2});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh{iS}(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh{iS}(:,:,ii) = mean(coh{1,1}.rel,3);
    end
end
mmCoh36 = cellfun(@(x) mean(x,3,'OmitNaN'), mCoh,'UniformOutput',0);
mstdCoh36 = cellfun(@(x) std(x,[],3,'OmitNaN'), mCoh,'UniformOutput',0);
mmRelCoh36 = cellfun(@(x) mean(x,3,'OmitNaN'), mRelCoh,'UniformOutput',0);
mstdRelCoh36 = cellfun(@(x) std(x,[],3,'OmitNaN'), mRelCoh,'UniformOutput',0);

%% 37
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\37\';
searchStr = {'base1','in';'base2','in';'alc1','in';'cloz1','in'};
% searchStr = {'base2','in';'alc2','in';'cloz4','in'};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1},searchStr{iS,2});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh{iS}(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh{iS}(:,:,ii) = mean(coh{1,1}.rel,3);
    end
end
mmCoh37 = cellfun(@(x) mean(x,3,'OmitNaN'), mCoh,'UniformOutput',0);
mstdCoh37 = cellfun(@(x) std(x,[],3,'OmitNaN'), mCoh,'UniformOutput',0);
mmRelCoh37 = cellfun(@(x) mean(x,3,'OmitNaN'), mRelCoh,'UniformOutput',0);
mstdRelCoh37 = cellfun(@(x) std(x,[],3,'OmitNaN'), mRelCoh,'UniformOutput',0);

%% Plot all coh for all conds
col = distinguishable_colors(6);
b = {'base1','base2','alc1','cloz1'};
for bi = 1:4
figure
hold on
c = 1;
for ii = 1:6
    shadedErrorBar(1:50,mmCoh36{1,bi}(ii,:),mstdCoh36{1,bi}(ii,:),{'--','Color',col(ii,:),'MarkerFaceColor',col(ii,:)},1)
    c = c+1;
end
title(['M36: Average Coherence across Animals ',b{bi}])
ylim([0 0.9])
figure
hold on
c = 1;
for ii = 1:6
    shadedErrorBar(1:50,mmCoh37{1,bi}(ii,:),mstdCoh37{1,bi}(ii,:),{'-','Color',col(ii,:),'MarkerFaceColor',col(ii,:)},1)
    c = c+1;
end
ylim([0 0.9])
title(['M37: Average Coherence across Animals ',b{bi}])
end
%% Plot cmbs 2,3,4,5 across conds
col = distinguishable_colors(6);
col(3,:) = [0 0.5 0];
col(6,:) = [1 0.5 0];
titles = {'SL-SR','SL-OL','SL-OR','SR-OL','SR-OR','OL-OR'};
for bi = 1:3%2:4
   f = figure;
   for fi = 1:6
       subtightplot(2,3,fi,[0.1,0.05],0.1,0.05)
       h = shadedErrorBar(1:2:100,mmCoh36{1,bi}(fi,:),mstdCoh36{1,bi}(fi,:),{'--','Color',col(fi,:),'MarkerFaceColor',col(fi,:)},1);
       h.edge(1,1).LineStyle = ':';
       h.edge(1,1).LineWidth = 1;
       h.edge(1,2).LineStyle = ':';
       h.edge(1,2).LineWidth = 1;
       hold on
       shadedErrorBar(1:2:100,mmCoh37{1,bi}(fi,:),mstdCoh37{1,bi}(fi,:),{'-','Color',col(fi,:),'MarkerFaceColor',col(fi,:)},1)
       ylim([0 1])
       title(titles{fi})
       if fi == 1 || fi == 4
           ylabel('Coherence')
       end
       if fi == 4 || fi == 5 || fi == 6
          xlabel('Frequency (Hz)') 
       end
   end
end
%% Collate data
[data36,samp36] = collateData('C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\36\',{'base2';'alc1';'cloz1'},{'coh'},'avg');
[data37,samp37] = collateData('C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\37\',{'base2';'alc1';'cloz1'},{'coh'},'avg');
%% Concatenate base1
base36 = cat(1,data36{1,1}{:});
alc36 = cat(1,data36{1,2}{:});
cloz36 = cat(1,data36{1,3}{:});
base37 = cat(1,data37{1,1}{:});
alc37 = cat(1,data37{1,2}{:});
cloz37 = cat(1,data37{1,3}{:});
%% Compare baselines
[h,p(1,:),pAdj(1,:)] = bulkT(base36,base37,0,'fdr');
[h,p(2,:),pAdj(2,:)] = bulkT(alc36,alc37,0,'fdr');
[h,p(3,:),pAdj(3,:)] = bulkT(cloz36,cloz37,0,'fdr');
%% Plot relCoh - one combo per panel; one condition per figure
cond = {'Baseline','Alcohol','Clozapine'};
for bi = 2:4
    figure
    for ii = 1:6
        subtightplot(2,3,ii,[0.1,0.05],0.1,0.05)
        b = barwitherr([mstdRelCoh36{1,bi}(ii,:);mstdRelCoh37{1,bi}(ii,:)]'.*100,[mmRelCoh36{1,bi}(ii,:);mmRelCoh37{1,bi}(ii,:)]'.*100);
        b(1).EdgeColor = col(ii,:);
        b(1).FaceColor = col(ii,:);
        b(2).EdgeColor = col(ii,:);
        b(2).FaceColor  = 'w';
        b(2).LineWidth = 2;
        set(gca,'XTickLabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})
        title(titles{ii})
        ylim([0 250])
        xlim([0.5 6.5])
        thisP = p(bi-1,1+(ii-1)*6:6*ii);
        pInd = [];
        pInd = logicFind(0.05,thisP,'<=');
        thisPadj = pAdj(bi-1,1+(ii-1)*6:6*ii);
        if ~isempty(pInd)
            for iP = 1:length(pInd)
                text(pInd(iP),200,'*')
            end
        end
        if ismember(ii,[1,4])
            ylabel('% of Average Coherence')
        end
        if ismember(ii,(4:6))
            xlabel('Frequency Band')
        end
    end
%     suptitle(['Relative Coherence: ',cond{bi-1}])
end
%% Plot relCoh - just looking at the trend level features across conditions
nums = [2,3;2,4;2,6];
titles = {'SL-OL: \alpha','SL-OL: \beta','SL-OL: h\gamma'};
pY = [1.3,1.3,1.1];
for ii = 1:size(nums,1)
    figure
    stds(:,:,ii) = [mstdRelCoh36{1,2}(nums(ii,1),nums(ii,2)),mstdRelCoh36{1,3}(nums(ii,1),nums(ii,2)),mstdRelCoh36{1,4}(nums(ii,1),nums(ii,2));mstdRelCoh37{1,2}(nums(ii,1),nums(ii,2)),mstdRelCoh37{1,3}(nums(ii,1),nums(ii,2)),mstdRelCoh37{1,4}(nums(ii,1),nums(ii,2))]';
    ms(:,:,ii) = [mmRelCoh36{1,2}(nums(ii,1),nums(ii,2)),mmRelCoh36{1,3}(nums(ii,1),nums(ii,2)),mmRelCoh36{1,4}(nums(ii,1),nums(ii,2));mmRelCoh37{1,2}(nums(ii,1),nums(ii,2)),mmRelCoh37{1,3}(nums(ii,1),nums(ii,2)),mmRelCoh37{1,4}(nums(ii,1),nums(ii,2))]';
    b = barwitherr(stds(:,:,ii),ms(:,:,ii));
    b(1).EdgeColor = col(nums(ii,1),:);
    b(1).FaceColor = col(nums(ii,1),:);
    b(2).EdgeColor = col(nums(ii,1),:);
    b(2).FaceColor  = 'w';
    b(2).LineWidth = 2;
    title(titles{ii})
    set(gca,'XTickLabel',{'Base','Alcohol','Clozapine'})
    ylabel('Normalized Coherence')
    text(1,pY(ii),'*','FontSize',20)
end 
%% Animals across conditions 
% Set up ylabels
ylab  = {'Base','Alc1','Alc2','Cloz1','Cloz2','Cloz3','Cloz4'};
%% 366
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\36\';
searchStr = {'366','in';};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1},searchStr{iS,2});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh366(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh366(:,:,ii) = mean(coh{1,1}.rel,3);
       sRelCoh366(:,:,ii) = std(coh{1,1}.rel,[],3);
    end
end
ord = [2,1,9,6,7,8,4,5];
% Plot all coherence
allSort(:,:,1) = squeeze(mCoh366(2,:,ord))';
figure
subplot(1,2,1)
imagesc(allSort(2:end,:,1))
set(gca,'YTickLabel',ylab)
title('36-6: Lesion')
colormap('viridis')
% Plot relative coherence
relSort(:,:,1) = squeeze(mRelCoh366(2,:,ord))';
figure
b = barwitherr(squeeze(sRelCoh366(2,:,ord)),squeeze(mRelCoh366(2,:,ord)));
set(gca,'XTickLabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})
xlabel('Frequency Band')
b(2).FaceColor = 'y'; b(3).FaceColor = 'y';
b(4).FaceColor = 'c'; b(5).FaceColor = 'c'; b(6).FaceColor = 'c'; b(7).FaceColor = 'c'; b(8).FaceColor = 'c';
%% 367
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\36\';
searchStr = {'367','in';};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1},searchStr{iS,2});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh367(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh367(:,:,ii) = mean(coh{1,1}.rel,3);
       sRelCoh367(:,:,ii) = std(coh{1,1}.rel,[],3);
    end
end
ord = [2,1,8,5,6,7,4,3];
% Plot all coherence
allSort(:,:,2) = squeeze(mCoh367(2,:,ord))';
figure
imagesc(allSort(2:end,:,2))
set(gca,'YTickLabel',ylab)
title('36-7: Lesion')
colormap('viridis')
% Plot relative coherence
relSort(:,:,2) = squeeze(mRelCoh367(2,:,ord))';
figure
b = barwitherr(squeeze(sRelCoh367(2,:,ord)),squeeze(mRelCoh367(2,:,ord)));
set(gca,'XTickLabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})
xlabel('Frequency Band')
b(2).FaceColor = 'y'; b(3).FaceColor = 'y';
b(4).FaceColor = 'c'; b(5).FaceColor = 'c'; b(6).FaceColor = 'c'; b(7).FaceColor = 'c'; b(8).FaceColor = 'c';
%% 374
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\37\';
searchStr = {'374','in';};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1},searchStr{iS,2});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh374(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh374(:,:,ii) = mean(coh{1,1}.rel,3);
       sRelCoh374(:,:,ii) = std(coh{1,1}.rel,[],3);
    end
end
ord = [1,2,6,7,8,4,5,3];
% Plot all coherence
allSort(:,:,3) = squeeze(mCoh374(2,:,ord))';
figure
imagesc(allSort(2:end,:,3))
set(gca,'YTickLabel',ylab)
title('37-4: Sham')
colormap('viridis')
% Plot relative coherence
relSort(:,:,3) = squeeze(mRelCoh374(2,:,ord))';
figure
b = barwitherr(squeeze(sRelCoh374(2,:,ord)),squeeze(mRelCoh374(2,:,ord)));
set(gca,'XTickLabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})
xlabel('Frequency Band')
b(2).FaceColor = 'y'; b(3).FaceColor = 'y';
b(4).FaceColor = 'c'; b(5).FaceColor = 'c'; b(6).FaceColor = 'c'; b(7).FaceColor = 'c'; b(8).FaceColor = 'c';
%% 375
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\37\';
searchStr = {'375','in';};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1},searchStr{iS,2});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh375(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh375(:,:,ii) = mean(coh{1,1}.rel,3);
       sRelCoh375(:,:,ii) = std(coh{1,1}.rel,[],3);
    end
end
ord = [1,2,6,7,8,4,5,3];
% Plot all coherence
allSort(:,:,4) = squeeze(mCoh375(2,:,ord))';
figure
imagesc(allSort(2:end,:,4))
set(gca,'YTickLabel',ylab)
title('37-5: Sham')
colormap('viridis')
% Plot relative coherence
relSort(:,:,4) = squeeze(mRelCoh375(2,:,ord))';
figure
b = barwitherr(squeeze(sRelCoh375(2,:,ord)),squeeze(mRelCoh375(2,:,ord)));
set(gca,'XTickLabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})
xlabel('Frequency Band')
b(2).FaceColor = 'y'; b(3).FaceColor = 'y';
b(4).FaceColor = 'c'; b(5).FaceColor = 'c'; b(6).FaceColor = 'c'; b(7).FaceColor = 'c'; b(8).FaceColor = 'c';
%% Plot average coherences
datMean(:,:,1) = mean(allSort(:,:,1:2),3);
datMean(:,:,2) = mean(allSort(:,:,3:4),3);
figure
imagesc(squeeze(datMean(2:end,:,1)))
set(gca,'YTickLabel',ylab,'CLim',[0.15 0.6],'XTickLabel',10:10:100)
title('Lesion Average')
colormap('viridis')
figure
imagesc(squeeze(datMean(2:end,:,2)))
set(gca,'YTickLabel',ylab,'CLim',[0.15 0.6],'XTickLabel',10:10:100)
title('Sham Average')
colormap('viridis')
clim([0.15 0.6])
%% Plot trending three variables across conditions
test1 = relMean(:,[3,4,6],1);
test2 = relMean(:,[3,4,6],2);
test3 = relStd(:,[3,4,6],1);
test4 = relStd(:,[3,4,6],2);
for ii = 1:3
    figure
    b = barwitherr([test3(2:end,ii),test4(2:end,ii)],[test1(2:end,ii),test2(2:end,ii)]);
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    title(titles{ii});
    set(gca,'XTickLabel',{'Base','Alc1','Alc2','Cloz1','Cloz2','Cloz3','Cloz4'})
    ylabel('Normalized Coherence')
end
%%
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\37\';
searchStr = {'base2';'alc1';'alc2';'cloz1';'cloz2';'cloz3';'cloz4'};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh37{iS}(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh37{iS}(:,ii) = squeeze(mean(coh{1,1}.rel(2,:,:),3));
       sRelCoh37{iS}(:,ii) = squeeze(std(coh{1,1}.rel(2,:,:),[],3));
    end
end
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\megan\processed\36\';
searchStr = {'base2';'alc1';'alc2';'cloz1';'cloz2';'cloz3';'cloz4'};
for iS = 1:size(searchStr,1)
    [files{iS}] = fileSearch(sdir,searchStr{iS,1});
    for ii = 1:size(files{iS},2)
       load(files{iS}{ii},'coh')
       mCoh36{iS}(:,:,ii) = coh{1,1}.mCxy;
       mRelCoh36{iS}(:,ii) = squeeze(mean(coh{1,1}.rel(2,:,:),3));
       sRelCoh36{iS}(:,ii) = squeeze(std(coh{1,1}.rel(2,:,:),[],3));
    end
end
%% Trend level vars all animals across conds
c = 1;
for k = [3,4,6]
    ml(c,:) = cellfun(@(x) mean(x(k,:),2,'omitNaN'),mRelCoh36,'UniformOutput',0);
    sl(c,:) = cellfun(@(x) std(x(k,:),[],2,'omitNaN'),mRelCoh36,'UniformOutput',0);
    ms(c,:) = cellfun(@(x) mean(x(k,:),2,'omitNaN'),mRelCoh37,'UniformOutput',0);
    ss(c,:) = cellfun(@(x) std(x(k,:),[],2,'omitNaN'),mRelCoh37,'UniformOutput',0);
    c = c+1;
end
%% Get p values for bunch of comparsions (L = lesion; S = sham; B = base; A = alcohol; C = clozapine)
c = 1;
for ii = [3,4,6]
    % LB-LA
    [~,p(1,c)] = ttest2(mRelCoh36{1}(ii,:),mRelCoh36{2}(ii,:));
    % LB-LC
    [~,p(2,c)] = ttest2(mRelCoh36{1}(ii,:),mRelCoh36{4}(ii,:));
    % LA-SB
    [~,p(3,c)] = ttest2(mRelCoh36{2}(ii,:),mRelCoh37{1}(ii,:));
    % LC-SB
    [~,p(4,c)] = ttest2(mRelCoh36{4}(ii,:),mRelCoh37{1}(ii,:));
    % SB-SA
    [~,p(5,c)] = ttest2(mRelCoh37{1}(ii,:),mRelCoh37{2}(ii,:));
    % SB-SC
    [~,p(6,c)] = ttest2(mRelCoh37{1}(ii,:),mRelCoh37{4}(ii,:));
    % LB-SB
    [~,p(7,c)] = ttest2(mRelCoh36{1}(ii,:),mRelCoh37{1}(ii,:));
    c = c+1;
end
grps = {[0.85 1.85],[.85 2.85],[1.85 1.15],[2.85 1.15],[1.15 2.15],[1.15 3.15],[.85 1.15]};
%%
titles = {'SL-OL: \alpha','SL-OL: \beta','SL-OL: h\gamma'};
for k = 1:3
   figure
   b = barwitherr([ss{k,[1,2,4]};sl{k,[1,2,4]}]'.*100,[ms{k,[1,2,4]};ml{k,[1,2,4]}]'.*100);
   b(1).FaceColor = 'w';
   b(2).FaceColor = 'k';
%    if k == 1
%        text(.9,130,'*','FontSize',20)
%        text(2.9,130,'*','FontSize',20)
%    end
%    if k == 2
%        text(.9,130,'*','FontSize',20)
%    end
%    if k == 3
%        text(.9,110,'*','FontSize',20)
%    end
   title(titles{k});
   set(gca,'XTickLabel',{'Base','Alc1','Cloz1'})
   ylabel('% of Average Coherence')
   thisP = mafdr(p(:,k),'BHFDR',1);
   pInd = logicFind(0.05,thisP,'<=');
   sigstar({grps{pInd}},[thisP(pInd)]);
end

%%
for ii = 1:7
    c = 1;
    for k = [3,4,6]
        [~,p(ii,c)] = ttest2(mRelCoh36{ii}(k,:),mRelCoh37{ii}(k,:));
        c=c+1;
    end
end
pAdj = mafdr(p(:,3),'BHFDR','true'); 
%% Plot average coherences
m36 = cellfun(@(x) mean(x(2,:,:),3,'omitnan'),mCoh36,'UniformOutput',0);
m37 = cellfun(@(x) mean(x(2,:,:),3,'omitnan'),mCoh37,'UniformOutput',0);
test1 = cat(1,m36{:});
test2 = cat(1,m37{:});
figure
imagesc(test1)
set(gca,'YTickLabel',ylab,'XTickLabel',10:10:100,'CLim',[0.2 0.6])
title('Lesion Average')
colormap('viridis')
figure
imagesc(test2)
set(gca,'YTickLabel',ylab,'XTickLabel',10:10:100,'CLim',[0.2 0.6])
title('Sham Average')
colormap('viridis')
%%
figure
h = shadedErrorBar(1:14,drink(1,1:14),drink(4,1:14),{'--k','Marker','o','MarkerSize',3},1);
h.edge(1,1).LineStyle = ':';
h.edge(1,2).LineStyle = ':';
hold on
shadedErrorBar(1:14,drink(2,1:14),drink(5,1:14),{'-k','Marker','o','MarkerFaceColor','k','MarkerSize',3},1)
h = shadedErrorBar(21:83,drink(1,16:end),drink(4,16:end),{'--k','Marker','o','MarkerSize',3},1);
h.edge(1,1).LineStyle = ':';
h.edge(1,2).LineStyle = ':';
shadedErrorBar(21:83,drink(2,16:end),drink(5,16:end),{'-k','Marker','o','MarkerFaceColor','k','MarkerSize',3},1);
xlabel('Days')
ylabel('Alcohol Consumed (gm/kgm)')
set(gca,'XTick',[0:7:14,27:7:83],'XTickLabel',{0:7:14,7:7:63})