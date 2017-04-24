% Examine first 3 minutes post stim relative to the average pre stim
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