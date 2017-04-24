function [T1,T2,fName,tsRest,tsBinge,nRest,nBinge] = collateData(sdir,searchStr,n,func)
%%
% Get file structure for each searchStr
[fileStruct] = fileSearch(sdir,searchStr);
% Go through each file structure and extract data from files therein
for fsi = 1:numel(fileStruct)
    disp(fsi)
    for fi = 1:size(fileStruct{fsi},1)
        % Splits filename at '_' and stores first part
        parts = strsplit(fileStruct{fsi}(fi).name,'_');
        fName{fsi,fi} = parts{1};
        load([sdir{1},fileStruct{fsi}(fi).name],'coh','psdTrls','relPower','powerCorrSort','LFPTs','trls','hist');
        % Stack overall PSD and PSD distributions
        %psdBinge{fsi}(:,:,fi) = psdTrls.event1.Overall;
        %psdRest{fsi}(:,:,fi) = psdTrls.event2.Overall;
        %relPSDbinge{fsi}(:,:,fi) = relPower.event1;
        %relPSDrest{fsi}(:,:,fi) = relPower.event2;
        %totalPSDbinge{fsi}(fi,:) = psdTrls.event1.totalPower;
        %totalPSDrest{fsi}(fi,:) = psdTrls.event2.totalPower;
        % Grab all power spectra
        powMatBinge{fsi,fi} = [];
        for ii = 1:length(psdTrls.event1.Pow)
            powMatBinge{fsi,fi} = cat(3,powMatBinge{fsi,fi},psdTrls.event1.Pow{1,ii});
        end
        powMatRest{fsi,fi} = [];
        for ii = 1:length(psdTrls.event2.Pow)
           powMatRest{fsi,fi} = cat(3,powMatRest{fsi,fi},psdTrls.event2.Pow{1,ii});
        end
        % Stack mean coherence
        meanCohBinge{fsi}(:,:,fi) = coh{1}.mCxy;
        %meanCohRest{fsi}(:,:,fi) = coh{2}.mCxy;
        cohMatBinge{fsi,fi} = coh{1}.Cxy;
        cohMatRest{fsi,fi} = coh{2}.Cxy;
        % Stack power correlations
        powCorrBinge{fsi,fi} = powerCorrSort{1}.masterCorr;
        powCorrRest{fsi,fi} = powerCorrSort{2}.masterCorr;
        % Stack number of windows
        nBinge{fsi}(fi) = size(psdTrls.event1.Pow,2);
        nRest{fsi}(fi) = size(psdTrls.event2.Pow,2);
        % Get sample start time-stamps
        tsBinge{fsi,fi} = trls{1,1}.sampleinfo(:,1);
        tsRest{fsi,fi} = trls{1,2}.sampleinfo(:,1);
        % On last iteration grab channel labels
        if fsi == numel(fileStruct) && fi == size(fileStruct{fsi},2)
            chanLabel = LFPTs.label;
            F = psdTrls.F;
            bands = hist.bands;
        end
        clear LFPTs trls relPower psdTrls coh stdPower hist powerEventComp powerCorrSort
        disp(fi)
    end
end
%% Normalize psd and coh by last n rest
% Get frequency band interval indices from indices in F
for j = 1:size(bands,1)
    bandInd(j,1) = find(F>=bands{j,2}(1),1);
    bandInd(j,2) = find(F<=bands{j,2}(2),1,'last');
end
for fsi = 1:numel(fileStruct)
    for fi = 1:size(fileStruct{fsi},1)
        if strcmpi(n,'all')
            inds = 1:nRest{fsi}(fi);
        end
        if strcmpi(n,'first')
            inds = 1:nBinge{fsi}(fi);
        end
        if strcmpi(n,'last')
            inds = nRest{fsi}(fi)-nBinge{fsi}(fi)+1:nRest{fsi}(fi);
        end
        if strcmpi(n,'rand')
            inds = randi(nRest{fsi}(fi),1,nBinge{fsi}(fi));
        end
        if strcmpi(func,'cat')
           for b = 1:size(bandInd,1)
                for c = 1:size(powMatBinge{fsi,fi},1)
                    relPsdRest{fsi,fi}(b,c,:) = sq(trapz(powMatRest{fsi,fi}(c,bandInd(b,1):bandInd(b,2),inds))./trapz(powMatRest{fsi,fi}(c,bandInd(1,1):bandInd(5,2),inds)));
                    relPsdBinge{fsi,fi}(b,c,:) = sq(trapz(powMatBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2),:))./trapz(powMatBinge{fsi,fi}(c,bandInd(1,1):bandInd(5,2),:)));
                end
                for c = 1:size(cohMatBinge{fsi,fi},1)
                    relCohRest{fsi,fi}(b,c,:) = sq(mean(cohMatRest{fsi,fi}(c,bandInd(b,1):bandInd(b,2),inds),2)./mean(cohMatRest{fsi,fi}(c,bandInd(1,1):bandInd(end,2),inds),2));
                    relCohBinge{fsi,fi}(b,c,:) = sq(mean(cohMatBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2),:),2)./mean(cohMatBinge{fsi,fi}(c,bandInd(1,1):bandInd(end,2),:),2));
                end
           end
           powCorrRestN{fsi,fi} = powCorrRest{fsi,fi}(inds,:);
        end
        if strcmpi(func,'mean')
           meanPsdRest{fsi,fi} = mean(powMatRest{fsi,fi}(:,:,inds),3);
           meanCohRest{fsi,fi} = mean(cohMatRest{fsi,fi}(:,:,inds),3);
           meanPowCorrRest{fsi,fi} = mean(powCorrRest{fsi,fi}(inds,:),1);
           meanPsdBinge{fsi,fi} = mean(powMatBinge{fsi,fi},3);
           meanCohBinge{fsi,fi} = mean(cohMatBinge{fsi,fi},3);
           meanPowCorrBinge{fsi,fi} = mean(powCorrBinge{fsi,fi},1);
           for b = 1:size(bandInd,1)
                for c = 1:size(meanPsdRest{fsi,fi},1)
                    relPsdRest{fsi}(b,c,fi) = trapz(meanPsdRest{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./trapz(meanPsdRest{fsi,fi}(c,bandInd(1,1):bandInd(5,2)));
                    relPsdBinge{fsi}(b,c,fi) = trapz(meanPsdBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./trapz(meanPsdBinge{fsi,fi}(c,bandInd(1,1):bandInd(5,2)));
                end
                for c = 1:size(meanCohRest{fsi,fi},1)
                    relCohRest{fsi}(b,c,fi) = mean(meanCohRest{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./mean(meanCohRest{fsi,fi}(c,bandInd(1,1):bandInd(end,2)));
                    relCohBinge{fsi}(b,c,fi) = mean(meanCohBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./mean(meanCohBinge{fsi,fi}(c,bandInd(1,1):bandInd(end,2)));
                end
            end
        end
    end
end

%% Build data array T with each column as a condition
% Assumes that if mean was used, there will be binge data and hard codes it
% in
if strcmpi(func,'mean')
    % Rest
    T1{1,1}(:,1) = [9.4;14.9;2.3;7.8;8.1;5.8;7.7;6.4;4.5;5.1;9;5.3;5.1;3;3.1;2.4;10;6.3;3.3;8;5;8.4;8.5;14.6];
    T1{1,2}(:,1) = [13.2;8.2;11.4;13.7;17.3;9.9;3.2;1.6;12.4;8.8;7.9;10.2];
    T1{1,3}(:,1) = [15.9;8.8;19.3;17.4;10.5;9.5;14.5;9.2;9.1];
    T1{1,4}(:,1) = [8.2;4.3;4;10.9;5.6;5.1;1.8;6.5];
    
    % Binge
    T2{1,1}(:,1) = [9.4;14.9;2.3;7.8;8.1;5.8;7.7;6.4;4.5;5.1;9;5.3;5.1;3;3.1;2.4;10;6.3;3.3;8;5;8.4;8.5;14.6];
    T2{1,2}(:,1) = [13.2;8.2;11.4;13.7;17.3;9.9;3.2;1.6;12.4;8.8;7.9;10.2];
    T2{1,3}(:,1) = [15.9;8.8;19.3;17.4;10.5;9.5;14.5;9.2;9.1];
    T2{1,4}(:,1) = [8.2;4.3;4;10.9;5.6;5.1;1.8;6.5];

    for fsi = 1:numel(fileStruct)
        for fi = 1:size(fileStruct{fsi},1)
            % Rest
            T1{1,fsi}(fi,2:21) = reshape(relPsdRest{fsi}(:,:,fi),[1,20]);
            T1{1,fsi}(fi,22:51) = reshape(relCohRest{fsi}(:,:,fi),[1,30]);
            T1{1,fsi}(fi,52:145) = reshape(meanPowCorrRest{fsi,fi},[1,94]);
            % Binge
            T2{1,fsi}(fi,2:21) = reshape(relPsdBinge{fsi}(:,:,fi),[1,20]);
            T2{1,fsi}(fi,22:51) = reshape(relCohBinge{fsi}(:,:,fi),[1,30]);
            T2{1,fsi}(fi,52:145) = reshape(meanPowCorrBinge{fsi,fi},[1,94]);
        end
    end
end
% Assumes if cat was used that this is binge vs not binge and inserts 1s
% and 0s accordingly
if strcmpi(func,'cat')
    for fsi = 1:numel(fileStruct)
        T1{1,fsi} = [];
        for fi = 1:size(fileStruct{fsi},1)
            chunkR = [zeros(size(relPsdRest{fsi,fi},3),1),reshape(relPsdRest{fsi,fi}(:,:,:),[20,size(relPsdRest{fsi,fi},3)])',reshape(relCohRest{fsi,fi}(:,:,:),[30,size(relCohRest{fsi,fi},3)])',powCorrRestN{fsi,fi}];
            chunkB = [ones(size(relPsdBinge{fsi,fi},3),1),reshape(relPsdBinge{fsi,fi}(:,:,:),[20,size(relPsdBinge{fsi,fi},3)])',reshape(relCohBinge{fsi,fi}(:,:,:),[30,size(relCohBinge{fsi,fi},3)])',powCorrBinge{fsi,fi}];
            % If all, create separate cells for each file
            if strcmpi(n,'all')
                T1{fsi,fi} = vertcat(chunkR,chunkB);
            else
                T1{1,fsi} = vertcat(T1{1,fsi},chunkR,chunkB);
            end
        end
    end
    T2 = [];
end
%% Subsets
% Get mean of baseline
% for c = 1:size(T1{1,1},2)
%     thisT1(:,c) = kMean(2,T1{1,1}(:,c));
%     thisT2(:,c) = kMean(2,T2{1,1}(:,c));
% end
% % Subset
% T1{2,1} = thisT1([1,2,3,4,5,6,11],:);
% T2{2,1} = thisT2([1,2,3,4,5,6,11],:);
% T1{2,2} = T1{1,2}([1,2,3,5,6,7,11],:);
% T2{2,2} = T2{1,2}([1,2,3,5,6,7,11],:);
% T1{2,3} = T1{1,3}([1,2,3,4,5,6,9],:);
% T2{2,3} = T2{1,3}([1,2,3,4,5,6,9],:);
% T1{2,4} = T1{1,4}([1,2,3,4,5,6,8],:);
% T2{2,4} = T2{1,4}([1,2,3,4,5,6,8],:);
% %% Do subtraction of binge from rest
% T3 = cellfun(@minus,T1,T2,'UniformOutput',false);
% T3{1,1}(:,1) = T1{1,1}(:,1);
%%
    % Subset
    % T1{2,1}(:,1) = kMean(2,T1{1,1}(:,1));
    % T1{2,2}(:,1) = [13.2;8.2;11.4;13.7;17.3;9.9;3.2;1.6;12.4;8.8;7.9];
    % T1{2,3}(:,1) = [15.9;8.8;19.3;17.4;10.5;9.5;14.5;9.2;9.1];
    % T1{3,4}(:,1) = [8.2;4.3;4;10.9;5.6;5.1;1.8;6.5];
    % Subset
    % T2{2,1}(:,1) = kMean(2,T2{1,1}(:,1));
    % T2{2,2}(:,1) = [13.2;8.2;11.4;13.7;17.3;9.9;3.2;1.6;12.4;8.8;7.9];
    % T2{2,3}(:,1) = [15.9;8.8;19.3;17.4;10.5;9.5;14.5;9.2;9.1];
    % T2{3,4}(:,1) = [8.2;4.3;4;10.9;5.6;5.1;1.8;6.5];
% if strcmpi(func,'mean')
%     for fsi = 1:numel(fileStruct)
%         for fi = 1:size(fileStruct{fsi},1)
%             % Takes all rests
%             if strcmpi(n,'all')
%                 meanPsdRestN{fsi,fi} = mean(powMatRest{fsi,fi}(:,:,:),3);
%                 meanCohRestN{fsi,fi} = mean(cohMatRest{fsi,fi}(:,:,:),3);
%                 meanPowCorrRestN{fsi,fi} = mean(powCorrRest{fsi,fi}(:,:),1);
%             end
%             % Uses last n rests equivalent to each nBinge
%             if strcmpi(n,'last')
%                 meanPsdRestN{fsi,fi} = mean(powMatRest{fsi,fi}(:,:,end-nBinge{fsi}(fi):end),3);
%                 meanCohRestN{fsi,fi} = mean(cohMatRest{fsi,fi}(:,:,end-nBinge{fsi}(fi):end),3);
%                 meanPowCorrRestN{fsi,fi} = mean(powCorrRest{fsi,fi}(end-nBinge{fsi}(fi):end,:),1);
%             end
%             % Uses randi to generate random row vector of integers
%             % between 1 and nRest of length nBinge
%             if strcmpi(n,'rand')
%                 randI = randi(nRest{fsi}(fi),1,nBinge{fsi}(fi));
%                 meanPsdRestN{fsi,fi} = mean(powMatRest{fsi,fi}(:,:,randI),3);
%                 meanCohRestN{fsi,fi} = mean(cohMatRest{fsi,fi}(:,:,randI),3);
%                 meanPowCorrRestN{fsi,fi} = mean(powCorrRest{fsi,fi}(randI,:),1);
%             end
%             meanPsdBinge{fsi,fi} = mean(powMatBinge{fsi,fi},3);
%             meanCohBinge{fsi,fi} = mean(cohMatBinge{fsi,fi},3);
%             meanPowCorrBinge{fsi,fi} = mean(powCorrBinge{fsi,fi},1);
%             for b = 1:size(bandInd,1)
%                 for c = 1:size(meanPsdRestN{fsi,fi},1)
%                     relPsdRestN{fsi}(b,c,fi) = trapz(meanPsdRestN{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./trapz(meanPsdRestN{fsi,fi}(c,bandInd(1,1):bandInd(5,2)));
%                     relPsdBinge{fsi}(b,c,fi) = trapz(meanPsdBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./trapz(meanPsdBinge{fsi,fi}(c,bandInd(1,1):bandInd(5,2)));
%                 end
%                 for c = 1:size(meanCohRestN{fsi,fi},1)
%                     relCohRestN{fsi}(b,c,fi) = mean(meanCohRestN{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./mean(meanCohRestN{fsi,fi}(c,bandInd(1,1):bandInd(end,2)));
%                     relCohBinge{fsi}(b,c,fi) = mean(meanCohBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2)))./mean(meanCohBinge{fsi,fi}(c,bandInd(1,1):bandInd(end,2)));
%                 end
%             end
%         end
%     end
% end
% if strcmpi(func,'cat')
%     for fsi = 1:numel(fileStruct)
%         for fi = 1:size(fileStruct{fsi},1)
%             % Takes all rests
%             if strcmpi(n,'all')
%                 inds = 1:size(nRest{fsi}(fi));
%             end
%             if strcmpi(n,'last')
%                 inds = 1:size(nBinge{fsi}(fi));
%             end
%             if strcmpi(n,'rand')
%                 inds = randi(nRest{fsi}(fi),1,nBinge{fsi}(fi));
%             end
%             for b = 1:size(bandInd,1)
%                 for c = 1:size(meanPsdRestN{fsi,fi},1)
%                     relPsdRest{fsi,fi}(c,b,:) = sq(trapz(powMatRest{fsi,fi}(c,bandInd(b,1):bandInd(b,2),inds))./trapz(powMatRest{fsi,fi}(c,bandInd(1,1):bandInd(5,2),inds)));
%                     relPsdBinge{fsi,fi}(c,b,:) = sq(trapz(powMatBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2),:))./trapz(powMatBinge{fsi,fi}(c,bandInd(1,1):bandInd(5,2),:)));
%                 end
%                 for c = 1:size(meanCohRestN{fsi,fi},1)
%                     relCohRest{fsi,fi}(c,b,:) = sq(mean(cohMatRest{fsi,fi}(c,bandInd(b,1):bandInd(b,2),inds),2)./mean(cohMatRest{fsi,fi}(c,bandInd(1,1):bandInd(end,2),inds),2));
%                     relCohBinge{fsi,fi}(c,b,:) = sq(mean(cohMatBinge{fsi,fi}(c,bandInd(b,1):bandInd(b,2),:),2)./mean(cohMatBinge{fsi,fi}(c,bandInd(1,1):bandInd(end,2),:),2));
%                 end
%             end
%         end
%     end
% end