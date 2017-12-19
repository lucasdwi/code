searchStr = {'5Hz_Pre','in';'5Hz_Post','in'};%'10Hz_Pre','in';'10Hz_Post','in'};
sdir = 'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\processed';
for sI = 1:size(searchStr,1)
    [files{sI}] = fileSearch(sdir,searchStr{sI,1},searchStr{sI,2});
    for fi = 1:size(files{sI},2)
       load(files{sI}{1,fi},'coh')
       cohMat{sI,fi}(:,:,:) = coh{1,1}.Cxy;
       bandMat{sI,fi}(:,:,:) = coh{1,1}.band;
    end
end
%% Normalize by total average
n = 0;
for rI = [1]
    for cI = 1:sum(~cellfun(@isempty,cohMat(rI,:)))
        if n ~= 0
            % Frequency average
            fM{rI,cI} = mean(cohMat{rI,cI}(:,:,end-n+1:end),3);
            fS{rI,cI} = std(cohMat{rI,cI}(:,:,end-n+1:end),[],3);
        else
            % Frequency average
            fM{rI,cI} = mean(cohMat{rI,cI},3);
            fS{rI,cI} = std(cohMat{rI,cI},[],3);
        end
        % All average
        aM{rI,cI} = mean(fM{rI,cI},2);
        aS{rI,cI} = std(fM{rI,cI},[],2);
        % Get number of samples
        s = size(cohMat{rI+1,cI},3);
        fN{rI,cI} = (cohMat{rI+1,cI}-repmat(fM{rI,cI},1,1,s))./repmat(fS{rI,cI},1,1,s);
        % Normalize by all
%         aN{rI,cI} = (cohMat{rI+1,cI}-repmat(aM{rI,cI},1,1,s))./repmat(aS{rI,cI},s,1000,1);
        
    end
end
%% Plot animal 6 band
for ii = 1
    for jj = 1:10
        mBand{jj} = mean(bandMat{ii,jj},3);
        sBand{jj} = std(bandMat{ii,jj},[],3);
        s = size(bandMat{ii+1,jj},3);
        zBand{jj} = (bandMat{ii+1,jj} - repmat(mBand{jj},1,1,s))./repmat(sBand{jj},1,1,s);
    end
end
%%
figure
for ii = 1:10
   subplot(5,2,ii)
   imagesc(squeeze(zBand{ii}(5,:,1:120)))
end
%%
test = 2*(1-normcdf(abs(squeeze(fN{9}(5,:,:))),0,1));
test(test>=0.05) = NaN;
figure
pcolor(test);
%%
for jj = 5
    figure
    for ii = 1:size(fN,2)
        subplot(5,2,ii)
%         test = 2*(1-normcdf(abs(squeeze(fN{ii}(jj,:,1:60))),0,1));
%         test(test>=0.05) = NaN;
%         pcolor(test)
            imagesc(squeeze(fN{ii}(jj,:,1:120)));
        colormap('viridis')
            caxis([-4 4])
    end
end
%% Test: perumatation test for significance
% First obtain p value from 'z' score calculated above
testP = 2*(1-normcdf(abs(fN{1,1}(1:60,1:500,1))));
figure
imagesc(testP')
% Randomize post stim data by frequency
for iterate = 1:1000
    randInd = randperm(500,500);
    randPost = cohMat{2,1}(:,randInd,1);
    randZ(:,:,iterate) = (randPost-repmat(fM{1,1}(:,1:500,1),360,1,1))./repmat(fS{1,1}(:,1:500,1),360,1,1);
end
%%
for ind1 = 1:size(randZ,1)
    for ind2 = 1:size(randZ,2)
        p(ind1,ind2) = sum(abs(squeeze(randZ(ind1,ind2,:)))>=abs(fN{1,1}(ind1,ind2,1)))/1000; 
    end
end
%%
test = fN{1,1}(:,1:500,1);
test(p>=0.05) = NaN;
%%
[h,critP,adjP] = fdr_bh(p,0.05,'dep');
%% Test plot
figure
subplot(2,2,1)
imagesc(301:360,1:0.2:200,cohMat{1,1}(301:360,1:500,1)')
subplot(2,2,3)
imagesc(1:60,1:0.2:200,cohMat{2,1}(1:60,1:500,1)')
subplot(2,2,2)
imagesc(1:60,1:0.2:200,fN{1,1}(1:60,1:500,1)')
subplot(2,2,4)
imagesc(1:60,1:0.2:200,aN{1,1}(1:60,1:500,1)')
%% Examine first 3 minutes post stim relative to the average pre stim
% 3 min = 180 sec = 36 5 sec trials
[~,fNames] = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\processed\'},{'pre','post'});
% Set number of trials from the beginning to look at in post file
n = 12;
% Takes advantage of automatic ordering from fileSearch to keep both fName
% arrays in same file order
for fi = 1:numel(fNames{1,1})
    load(fNames{1,1}{1,fi},'psdTrls','coh')
    % Grab mean and std from pre
    preMeanPow = psdTrls.event1.Overall;
    preSdPow = psdTrls.event1.OverallStd;
    preMeanCoh = coh{1,1}.mCxy;
    preSdCoh = coh{1,1}.sdCxy;
    load(fNames{1,2}{1,fi},'psdTrls','coh')
    % Grab first n psds from post
    postPow = cat(3,psdTrls.event1.Pow{1,1:n});
    postCoh = cat(3,coh{1,1}.Cxy(:,:,1:n));
    % Normalize each post trial by mean and std of pre
    for ni = 1:n
        normPostPow(:,:,ni,fi) = (postPow(:,:,ni)-preMeanPow)./preSdPow;
        normPostCoh(:,:,ni,fi) = (postCoh(:,:,ni)-preMeanCoh)./preSdCoh;
    end
    % end
    % % Get average 'z' scores
    % normPostPowAvg = mean(normPostPow,4);
    % normPostCohAvg = mean(normPostCoh,4);
    % % Get standard deviations of 'z' scores
    % normPostPowSd = std(normPostPow,[],4);
    % normPostCohSd = std(normPostCoh,[],4);
    % Get p-values from 'z'-scores
    normPostPowP = 2*(1-normcdf(abs(normPostPow)));
    normPostCohP = 2*(1-normcdf(abs(normPostCoh)));
    % normPostPowP = 2*(1-normcdf(abs(normPostPowAvg)));
    % normPostCohP = 2*(1-normcdf(abs(normPostCohAvg)));
    % Apply BH Bonferroni correction
    [~,~,pAdjPow] = fdr_bh(normPostPowP,0.05,'pdep');
    [~,~,pAdjCoh] = fdr_bh(normPostCohP,0.05,'pdep');
    % Set all pAdjSign > 0.05 to zero
    pAdjPow(pAdjPow > 0.05) = 0;
    pAdjCoh(pAdjCoh > 0.05) = 0;
    % Apply signs from data to pAdj to indicate direction
    pAdjSignPow = pAdjPow.*sign(normPostPow);
    pAdjSignCoh = pAdjCoh.*sign(normPostCoh);
end
%%
for fi = 1:numel(fNames{1,1})
    % Plot Power
    cMaxPow = max(max(max(max(pAdjSignPow))));
    cMinPow = min(min(min(min(pAdjSignPow))));
    figure
    for si = 1:4
        subplot(2,2,si)
        imagesc(sq(pAdjSignPow(si,:,:,fi)),[cMinPow cMaxPow])
    end
    % Plot Coherence
    cMaxCoh = max(max(max(max(pAdjSignCoh))));
    cMinCoh = min(min(min(min(pAdjSignCoh))));
    figure
    for si = 1:size(pAdjSignCoh,1)
        subplot(3,2,si)
        imagesc(sq(pAdjSignCoh(si,:,:,fi)),[cMinCoh cMaxCoh])
    end
end