% [~,fNames] = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\'},{'base'});
[~,fNames] = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\'},{'binge'});
% [~,fNames] = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\'},{'H10','H13','H14','H15','I1a','I2','I3','I4','I6','I8','I11','I12'});
for fsi = 1:size(fNames,2)
    disp(num2str(fsi))
    for fi = 1:size(fNames{1,fsi},2)
        disp(num2str(fi))
        load(fNames{1,fsi}{1,fi},'psdTrls','hist','coh')
        nBinge(fsi,fi) = size(psdTrls.event1.Pow,2);
        nRest(fsi,fi) = size(psdTrls.event2.Pow,2);
        for j = 1:size(hist.bands,1)
            bandInd(j,1) = find(psdTrls.F>=hist.bands{j,2}(1),1);
            bandInd(j,2) = find(psdTrls.F<=hist.bands{j,2}(2),1,'last');
        end
        % Grab all power and normalize
        for ii = 1:size(psdTrls.event1.Pow,2)
            for c = 1:4
                bingeTot(ii,c) = trapz(psdTrls.event1.Pow{1,ii}(c,1:45));
            end
            bingeRel(:,:,ii) = psdTrls.event1.Pow{2,ii}./repmat(bingeTot(ii,:),6,1);
            bingePow(ii,:) = reshape(bingeRel(:,:,ii),1,24);
        end
        for ii = 1:size(psdTrls.event2.Pow,2)
            for c = 1:4
                restTot(ii,c) = trapz(psdTrls.event2.Pow{1,ii}(c,1:45));
            end
            restRel(:,:,ii) = psdTrls.event2.Pow{2,ii}./repmat(restTot(ii,:),6,1);
            restPow(ii,:) = reshape(restRel(:,:,ii),1,24);
        end
        % Grab all coherence and normalize
        % Band x cmb x trial
        for bi = 1:size(bandInd,1)
            bingeAvgCoh(bi,:,:) = mean(coh{1,1}.Cxy(:,bandInd(bi,1):bandInd(bi,2),:),2);
            restAvgCoh(bi,:,:) = mean(coh{1,2}.Cxy(:,bandInd(bi,1):bandInd(bi,2),:),2);
        end
        % cmb x trial
        allBingeCoh = squeeze(mean(coh{1,1}.Cxy(:,1:45,:),2));
        allRestCoh = squeeze(mean(coh{1,2}.Cxy(:,1:45,:),2));
        for ii = 1:size(allBingeCoh,2)
            bingeRelCoh(:,:,ii) = bingeAvgCoh(:,:,ii)./repmat(allBingeCoh(:,ii)',6,1);
            bingeCoh(ii,:) = reshape(bingeRelCoh(:,:,ii),1,36);
        end
        for ii = 1:size(allRestCoh,2)
            restRelCoh(:,:,ii) = restAvgCoh(:,:,ii)./repmat(allRestCoh(:,ii)',6,1);
            restCoh(ii,:) = reshape(restRelCoh(:,:,ii),1,36);
        end
        % Combine power and coherence
        allData{fsi,fi} = [bingePow,bingeCoh;restPow,restCoh];
        resp{fsi,fi} = [ones(size(bingePow,1),1);zeros(size(restPow,1),1)];
        clearvars -except fNames fi allData resp nBinge nRest fsi
    end
end