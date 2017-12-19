%%
sdir = 'C:\Users\Pythia\Documents\GreenLab\data\emily\csv\EK\';
postFiles = fileSearch(sdir,'post','in')';
for ii = 1:size(postFiles,1)
    % Extracting and concatenating columns B and H
    postData{ii,1} = [xlsread(postFiles{ii,1},'B:B'),[NaN;xlsread(postFiles{ii,1},'H:H')]];
end
% Find shortest file
[r,~] = cellfun(@size,postData);
minPost = min(r);
%
preFiles = fileSearch(sdir,'pre','in')';
for ii = 1:size(preFiles,1)
    % Extracting and concatenating columns B and H
    preData{ii,1} = [xlsread(preFiles{ii,1},'B:B'),[NaN;xlsread(preFiles{ii,1},'H:H')]];
end
% Find shortest file
[r,~] = cellfun(@size,preData);
minPre = min(r);
for ii = 1:size(postFiles,1)
    postDistance(ii,1) = sum(postData{ii,1}(1:minPost,2),'omit');
    preDistance(ii,1) = sum(preData{ii,1}(1:minPre,2),'omit');
end
%%
preDistance = reshape(preDistance,6,10);
postDistance = reshape(postDistance,6,10);
%%
for ii = 1:size(postData,1)
   for jj = 1:6
      postBinSum(ii,jj) = sum(postData{ii}(18001*jj-18000:18001*jj,2),'omit'); 
   end
%    postBinSum(ii,jj+1) = sum(postData{ii}(18001*(jj+1)-18000:end,2),'omit');
end
%%
postBinSum = reshape(postBinSum',6,10,6);
%%
figure
for ii = 1:6
    subplot(2,3,ii)
    scatterErr(1:6,mean(postBinSum(:,group==1,ii),2),std(postBinSum(:,group==1,ii),[],2),0)
    scatterErr(1.1:6.1,mean(postBinSum(:,group==0,ii),2),std(postBinSum(:,group==0,ii),[],2),0,[0.8 0.8 0.8])
    title(['Day ',num2str(ii)])
    xlim([0.5 5.5])
    ylim([0 5000])
    ylabel('Distance')
    xlabel('10 Minute Bin')
end
%%
h(1) = scatterErr(1:size(preDistance,1),mean(preDistance(:,group==1),2),std(preDistance(:,group==1),[],2),1,'k');
h(2) = scatterErr(1.1:size(preDistance,1)+0.1,mean(preDistance(:,group==0),2),std(preDistance(:,group==0),[],2),0,[0.8 0.8 0.8]);
title('Pre-Injection')
set(gca,'xtick',1:size(preDistance,1))
xlabel('Day')
ylabel('Distance (cm)')
legend([h(1) h(2)],{'Lesion','Sham'},'Location','southeast')

h(1) = scatterErr(1:size(postDistance,1),mean(postDistance(:,group==1),2),std(postDistance(:,group==1),[],2),1,'k');
h(2) = scatterErr(1.1:size(postDistance,1)+0.1,mean(postDistance(:,group==0),2),std(postDistance(:,group==0),[],2),0,[0.8 0.8 0.8]);
title('Post-Injection')
set(gca,'xtick',1:size(postDistance,1))
xlabel('Day')
ylabel('Distance (cm)')
legend([h(1) h(2)],{'Lesion','Sham'},'Location','southeast')
%%
endPre = xlsread('C:\Users\Pythia\Documents\GreenLab\data\emily\startStopTimes.csv','C:C');
startPost = xlsread('C:\Users\Pythia\Documents\GreenLab\data\emily\startStopTimes.csv','F:F');
sdir = 'C:\Users\Pythia\Documents\GreenLab\data\emily\csv\';
files = fileSearch(sdir,'.csv','in')';
data = cell(size(files,1),1);
indexEndPre = zeros(size(files,1),1);
indexStartPost = zeros(size(files,1),1);
for ii = 1:size(files,1)
    % Prep filename
    filename = files{ii,1};
    % Extracting and concatenating columns B and H
    data{ii,1} = [xlsread(filename,'B:B'),xlsread(filename,'H:H')];
    % Find end of pre and start of post indices
    indexEndPre(ii,1) = nearest_idx3(endPre(ii,1),data{ii,1}(:,1));
    indexStartPost(ii,1) = nearest_idx3(startPost(ii,1),data{ii,1}(:,1));
end
%% Find shortest file
[r,~] = cellfun(@size,data);
[~,shortInd] = min(r);
shortestTime = data{shortInd,1}(end,1)-startPost(shortInd,1);
%% Find distance traveled pre and post injection
preDistance = zeros(size(files,1),1);
postDistance = zeros(size(files,1),1);
for ii = 1:size(files,1)
    preDistance(ii,1) = sum(data{ii,1}(1:indexEndPre,2),'omit');
    % Apply shortestTime to each startPost
    thisEndPost = startPost(ii,1) + shortestTime;
    postDistance(ii,1) = sum(data{ii,1}(indexStartPost(ii,1):...
        nearest_idx3(thisEndPost,data{ii,1}(:,1)),2),'omit');
end
%%
% load('processed.mat')
pre = reshape(preDistance,4,15);
h(1) = scatterErr(1:4,mean(pre(:,group==1),2),std(pre(:,group==1),[],2),1,'k');
h(2) = scatterErr(1.1:4.1,mean(pre(:,group==0),2),std(pre(:,group==0),[],2),0,[0.8 0.8 0.8]);
title('Baseline')
set(gca,'xtick',1:4)
xlabel('Day')
ylabel('Distance (cm)')
legend([h(1) h(2)],{'THC','Vehicle'})
%%
post = reshape(postDistance,4,15);
h(1) = scatterErr(1:4,mean(post(:,group==1),2),std(post(:,group==1),[],2),1,'k');
h(2) = scatterErr(1.1:4.1,mean(post(:,group==0),2),std(post(:,group==0),[],2),0,[0.8 0.8 0.8]);
set(gca,'xtick',1:4)
xlabel('Day')
ylabel('Distance (cm)')
title('Post')
legend([h(1) h(2)],{'THC','Vehicle'})

