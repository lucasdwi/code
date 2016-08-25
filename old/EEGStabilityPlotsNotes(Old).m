% Days between recordings of H10, H13, H14, H15, I11, I12, I1, I2, I3, I4, I6, I8
day = [18,46,16,12,1,13,14,34,2,71,67,36];
first = [1:2:size(T.Base,1)];
second = [2:2:size(T.Base,1)];
for s = 1:length(first)
    sub(s,:) = table2array(T.Base(second(s),8:end)) - table2array(T.Base(first(s),8:end));
    absSub(s,:) = abs(table2array(T.Base(second(s),8:end)) - table2array(T.Base(first(s),8:end)));
end
%% Get NACL theta and NASR-NACR low gamma
for s = 1:length(first)
        theta_sub(s) = table2array(T.Base(second(s),18)) - table2array(T.Base(first(s),18));
        lgamma_sub(s) = table2array(T.Base(second(s),51)) - table2array(T.Base(first(s),51));
end
%% Pull out power for each band
% Initialize absolute power matrix
% absM= [];
% for ii = 1:5
%     % Get column indices
%     b{ii} = [ii:5:20];
%     % Initialize and build matrix of days and powers
%     bm{ii} = [];
%     for j = 1:4
%         bm{ii} = vertcat(bm{ii},[sub(:,b{ii}(j)),day']);
%     end
%     % Get mean abosulte difference across all channels
%     for j = 1:size(absSub,1)
%         absM(j,ii) = mean(absSub(j,[b{ii}]));
%         bavg(j,ii) = mean(sub(j,[b{ii}]));
%     end
% end
%% Pull out coherence for each band
% absCoh = [];
% for ii = 1:5
%     cb{ii} = [20+ii:5:50];
%     % Get mean absolute difference across all channel pairs
%     for j = 1:size(absSub,1)
%         absCoh(j,ii) = mean(absSub(j,[cb{ii}]));
%     end
% end
%% Coherence difference per channel pair
figure
mark = {'+','s','o','*','d'};
pairNames = {'NASL-NASR','NASL-NACL','NASL-NACR','NASR-NACL','NASR-NACR','NACL-NACR'};
for ii = 1:6
    subplot(2,3,ii)
    % Get indices for five columns of each channel pair
    pairInd{ii} = [21+5*(ii-1):25+5*(ii-1)];
    % Raw difference
    thisP = sub(:,pairInd{ii});
    % Absoulate difference
    %thisP = abs(sub(:,pairInd{ii}));
    % Run fitlm; band X pair
    hold on
    for m = 1:5
        md4{m,ii} = fitlm(day,thisP(:,m));
    end
    hold on
    % Raw difference
    plot(day,sub(:,pairInd{ii}),'o')
    % Absolute differnce
    %plot(day,abs(sub(:,pairInd{ii})),'o')
    h = lsline;
    xlabel('Number of Days'); ylabel('Change in Coherence Index');
    title(['Change in Coherence Index in ',pairNames{ii}]);
    legend({strcat('Theta p =',num2str(md4{1,ii}.Coefficients.pValue(2),2),'; r2 =',num2str(md4{1,ii}.Rsquared.Ordinary,2))...
    strcat('Alpha p =',num2str(md4{2,ii}.Coefficients.pValue(2),2),'; r2 =',num2str(md4{2,ii}.Rsquared.Ordinary,2))...
    strcat('Beta p =',num2str(md4{3,ii}.Coefficients.pValue(2),2),'; r2 =',num2str(md4{3,ii}.Rsquared.Ordinary,2))...
    strcat('L Gamma p =',num2str(md4{4,ii}.Coefficients.pValue(2),2),'; r2 =',num2str(md4{4,ii}.Rsquared.Ordinary,2))...
    strcat('H Gamma p =',num2str(md4{5,ii}.Coefficients.pValue(2),2),'; r2 =',num2str(md4{5,ii}.Rsquared.Ordinary,2))...
    },'Location','north');
end
%% Coherence absolute difference plot average across channel pairs
c = {[0 0.447 0.7410];[0.85 0.325 0.098];[0.929 0.694 0.125];[0.494 0.184 0.556];[0.4660 0.6740 0.1880]};
figure; hold on
for ii = 1:5
    plot(day,absCoh(:,ii),'o','Color',c{ii});
    md3{ii} = fitlm(day,absCoh(:,ii));
end
h = lsline;
xlabel('Number of Days'); ylabel('Change in Coherence Index'); 
title('Change in coherence index over time');
legend({strcat('Theta p =',num2str(md3{1,1}.Coefficients.pValue(2),2),'; r2 =',num2str(md3{1,1}.Rsquared.Ordinary,2))...
    strcat('Alpha p =',num2str(md3{1,2}.Coefficients.pValue(2),2),'; r2 =',num2str(md3{1,2}.Rsquared.Ordinary,2))...
    strcat('Beta p =',num2str(md3{1,3}.Coefficients.pValue(2),2),'; r2 =',num2str(md3{1,3}.Rsquared.Ordinary,2))...
    strcat('L Gamma p =',num2str(md3{1,4}.Coefficients.pValue(2),2),'; r2 =',num2str(md3{1,4}.Rsquared.Ordinary,2))...
    strcat('H Gamma p =',num2str(md3{1,5}.Coefficients.pValue(2),2),'; r2 =',num2str(md3{1,5}.Rsquared.Ordinary,2))...
    },'Location','north');
%% Raw difference plot averaged across channels
c = {[0 0.447 0.7410];[0.85 0.325 0.098];[0.929 0.694 0.125];[0.494 0.184 0.556];[0.4660 0.6740 0.1880]};
figure; hold on
for ii = 1:5
%     plot(bm{ii}(:,2),bm{ii}(:,1).*100,'o','Color',c{ii});
%     md1{ii} = fitlm(bm{ii}(:,2),bm{ii}(:,1).*100);
    plot(day,bavg(:,ii).*100,'o','Color',c{ii});
    md1{ii} = fitlm(day,bavg(:,ii).*100);
end
h = lsline;
xlabel('Number of Days'); ylabel('Average Percent Difference'); 
title('Change in percent power over time');
legend({strcat('Theta p =',num2str(md1{1,1}.Coefficients.pValue(2),2),'; r2 =',num2str(md1{1,1}.Rsquared.Ordinary,2))...
    strcat('Alpha p =',num2str(md1{1,2}.Coefficients.pValue(2),2),'; r2 =',num2str(md1{1,2}.Rsquared.Ordinary,2))...
    strcat('Beta p =',num2str(md1{1,3}.Coefficients.pValue(2),2),'; r2 =',num2str(md1{1,3}.Rsquared.Ordinary,2))...
    strcat('L Gamma p =',num2str(md1{1,4}.Coefficients.pValue(2),2),'; r2 =',num2str(md1{1,4}.Rsquared.Ordinary,2))...
    strcat('H Gamma p =',num2str(md1{1,5}.Coefficients.pValue(2),2),'; r2 =',num2str(md1{1,5}.Rsquared.Ordinary,2))...
    },'Location','north');
%% Absolute difference plot averaged across channels
c = {[0 0.447 0.7410];[0.85 0.325 0.098];[0.929 0.694 0.125];[0.494 0.184 0.556];[0.4660 0.6740 0.1880]};
figure; hold on
for ii = 1:5
    plot(day,absM(:,ii).*100,'o','Color',c{ii});
    md2{ii} = fitlm(day,absM(:,ii).*100);
end
h = lsline;
xlabel('Number of Days'); ylabel('Average Percent Absolute Difference'); 
title('Absolute change in percent power over time');
legend({strcat('Theta p =',num2str(md2{1,1}.Coefficients.pValue(2),2),'; r2 =',num2str(md2{1,1}.Rsquared.Ordinary,2))...
    strcat('Alpha p =',num2str(md2{1,2}.Coefficients.pValue(2),2),'; r2 =',num2str(md2{1,2}.Rsquared.Ordinary,2))...
    strcat('Beta p =',num2str(md2{1,3}.Coefficients.pValue(2),2),'; r2 =',num2str(md2{1,3}.Rsquared.Ordinary,2))...
    strcat('L Gamma p =',num2str(md2{1,4}.Coefficients.pValue(2),2),'; r2 =',num2str(md2{1,4}.Rsquared.Ordinary,2))...
    strcat('H Gamma p =',num2str(md2{1,5}.Coefficients.pValue(2),2),'; r2 =',num2str(md2{1,5}.Rsquared.Ordinary,2))...
    },'Location','north');