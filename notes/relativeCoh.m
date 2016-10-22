files = dir('*.mat*');
%
bands = {'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
cmb = nchoosek(1:4,2);
for ii = 1:numel(files)
    load(files(ii).name)
    for j = 1:size(bands,1)
        % Create frequency band intervals from indices in freq
        bandInd(j,1) = find(coh.freq>=bands{j,2}(1),1);
        bandInd(j,2) = find(coh.freq<=bands{j,2}(2),1,'last');
        for k = 1:size(cmb,1)
            bandCoh(j,k) = mean(coh.mCxy(k,bandInd(j,1):bandInd(j,2),:));
            relCoh(j,k) = bandCoh(j,k)./mean(coh.mCxy(k,:,:));
        end
    end
    coh.rel = relCoh;
    save(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisBingeTest\',files(ii).name(1:end-4)],'LFPTs','trls','relPower','psdTrls','coh','hist','powerCorrSort')
    clearvars -except files ii bands cmb
end