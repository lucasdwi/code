for ii = 1:12
    c = 1;
    for k = 1+((ii-1)*20):ii*20
        load([num2str(k),'.mat'],'x','y','a','predY')
        if size(x,1) == 2
           x = (0:1/size(predY,1):1);
           y = x;
        end
        xAll{ii}(c,:) = x;
        yAll{ii}(c,:) = y;
        aAll(ii,c) = a;
        c = c+1;
    end
end
%%
mA = mean(aAll,2);
sA = std(aAll,[],2);
mX = cellfun(@(x) mean(x,1),xAll,'UniformOutput',0);
sX = cellfun(@(x) std(x,[],1),xAll,'UniformOutput',0);
mY = cellfun(@(y) mean(y,1),yAll,'UniformOutput',0);
sY = cellfun(@(y) std(y,[],1),yAll,'UniformOutput',0);
z = (0.5-mA)./sA;
%%
titles = {'H10','H13','H14','H15','I11','I12','I1','I2','I3','I4','I6','I8'};
figure
for k = 1:12
    subplot(3,4,k)
    hold on
    for ii = 1:20
        plot(xAll{1,k}(ii,:),yAll{1,k}(ii,:),'b','LineWidth',0.1)
    end
    plot(mX{1,k},mY{1,k},'r','LineWidth',2)
    text(0.5,0.1,['z = ',num2str(abs(round(z(k),2)))]);
    title(num2str(k))
end
%%
for ii = 1:12
    c = 1;
    for k = 1+((ii-1)*20):ii*20
        load([num2str(k),'.mat'],'allX','allY','allA','hist')
%         for jj = 1:42
%             if size(allX{jj},1) == 2
%                 allX{jj} = (0:1/size(hist.testInd,1):1)';
%                 allY{jj} = allX{jj};
%             end
%         end
        xAll{ii}(c,:) = allX;
        yAll{ii}(c,:) = allY;
        aAll{ii}(c,:) = allA;
        c = c+1;
    end
end
%% Go through xAll and yAll, look for 'empty' models, and replace with correct dimensions
for ii = 1:12
   thisSize = cellfun(@(x) size(x,1),xAll{1,ii},'UniformOutput',0);
   for k = 1:42
        u = unique(cell2mat(thisSize(:,k)));
        if size(u,1) >= 2
            maxS = max(u);
            inds = logicFind(2,cell2mat(thisSize(:,k)),'==');
            xAll{1,ii}(inds,k) = {(0:1/(maxS-1):1)'};
            yAll{1,ii}(inds,k) = {(0:1/(maxS-1):1)'};
        end
   end
end
%%
mA = [];
mA = cellfun(@(x) mean(x,1),aAll,'UniformOutput',0);
for ii = 1:12
   test(ii,:) = mA{1,ii}; 
end
figure
imagesc(test); colormap('viridis')
%% Baseline - self
for ii = 1:12
    testx{ii} = cat(2,xAll{1,ii}{:,ii});
    testy{ii} = cat(2,yAll{1,ii}{:,ii});
end
titles = {'H10','H13','H14','H15','I11','I12','I1','I2','I3','I4','I6','I8'};
figure
for ii = 1:12
   subplot(3,4,ii)
   for k = 1:20
       hold on
       plot(testx{ii}(:,k),testy{ii}(:,k),'b','LineWidth',0.1)
   end
   plot(mean(testx{ii},2),mean(testy{ii},2),'r','LineWidth',2)
   title(titles{ii})
end
%% Dep 24 - self
for ii = 1:12
    testx{ii} = cat(2,xAll{1,ii}{:,ii+12});
    testy{ii} = cat(2,yAll{1,ii}{:,ii+12});
end
titles = {'H10','H13','H14','H15','I11','I12','I1','I2','I3','I4','I6','I8'};
figure
for ii = 1:12
   subplot(3,4,ii)
   for k = 1:20
       hold on
       plot(testx{ii}(:,k),testy{ii}(:,k),'b','LineWidth',0.1)
   end
   plot(mean(testx{ii},2),mean(testy{ii},2),'r','LineWidth',2)
   title(title{ii})
end
%% Dep 48 - self
dInds = 25:33;
c = 1;
for ii = [1,2,3,5,6,7,9,10,11]
    testx{c} = cat(2,xAll{1,ii}{:,dInds(c)});
    testy{c} = cat(2,yAll{1,ii}{:,dInds(c)});
    c = c+1;
end
titles = {'H10','H13','H14','H15','I11','I12','I1','I2','I3','I4','I6','I8'};
figure
c = 1;
for ii = [1,2,3,5,6,7,9,10,11]
   subplot(3,3,c)
   for k = 1:20
       hold on
       plot(testx{c}(:,k),testy{c}(:,k),'b','LineWidth',0.1)
   end
   plot(mean(testx{c},2),mean(testy{c},2),'r','LineWidth',2)
   title(titles{ii})
   c = c+1;
end
%% Chow - self
dInds = 34:42;
c = 1;
for ii = [1,2,3,5,6,7,8,10,11]
   chowX{c} = cat(2,xAll{1,ii}{:,dInds(c)});
   chowY{c} = cat(2,yAll{1,ii}{:,dInds(c)});
   c = c+1; 
end
titles = {'H10','H13','H14','H15','I11','I12','I1','I2','I3','I4','I6','I8'};
figure
c = 1;
for ii = [1,2,3,5,6,7,8,10,11]
   subplot(3,3,c)
   for k = 1:20
       hold on
       plot(chowX{c}(:,k),chowY{c}(:,k),'b','LineWidth',0.1)
   end
   plot(mean(chowX{c},2),mean(chowY{c},2),'r','LineWidth',2)
   title(titles{ii})
   c = c+1;
end