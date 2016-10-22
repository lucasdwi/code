function [STDCorr,MeanCorr,TWCorr,powerCorrSort,freqRange,notchInd,numCmb,varCmb] = powerCorr(powCorrTFR,bands,cfg)
%% Use on ft_freqanalysis output to sort frequencies per trial per channel
%Uses powerspectrum data files from TFR to compare contiguous frequency
%bands across channels based on comparison iteration. Current iterations
%with respective comparisons are listed as follows:
% e.g. Iteration = Frequency A x Frequency B (includes self comparisons)
% Iteration 1 = theta x alpha
% Iteration 2 = alpha x beta
% Iteration 3 = beta x low gamma
% Iteration 4 = low gamma x high gamma

% Due to utilization of nchoosek, self-comparisons within bands are
% accounted for, necessitating for only 4 iterations between the 5
% frequency bands. For additional frequencies or comparisions, 
% iterations=frequencybands-1.

%% Initialize frequency ranges and notch range
for iteration=1:length(bands)-1;
    % Tune for specific iterations if desired
    % iteration=1;
    
    %Resets data struct for each new data set
    if iteration == 1      
        powerCorrSort = [];
    end
    % Internal iteration cycle input, unnecessary when using function
    % iteration=4;
    % Setting indices determined by iteration
    ind1 = iteration;
    ind2 = ind1+1;
    % Setting ceiling for iteration indices
    if ind1 == size(bands,1);
        ind2 = ind1;
        display('Error: Iteration number equals frequency bands and will cause erroenous autocorrelations')
    end
    % Setting indices for Frequency A based on iteration and bands input
    bandInd1 = find(powCorrTFR{1,ind1}.freq>=bands{ind1,2}(1),1);
    bandInd2 = find(powCorrTFR{1,ind1}.freq<=bands{ind1,2}(2),1,'last');
    % Setting indices for Frequency B based on iteration and bands input
    bandInd3 = find(powCorrTFR{1,ind1}.freq>=bands{ind2,2}(1),1);
    bandInd4 = find(powCorrTFR{1,ind1}.freq<=bands{ind2,2}(2),1,'last');
    
    % Setting indices for total frequency array
    freqstart = find(powCorrTFR{1,ind1}.freq>=bands{1,2}(1),1);
    freqend = find(powCorrTFR{1,ind1}.freq<=bands{end,2}(2),1,'last');
    % Create indices for total frequencies between first bands to last bands,
    % set at function level as 4Hz:90Hz
    freqRange = freqstart:freqend;
    F = powCorrTFR{1,ind1}.freq;
    %Create notch index
    notchInd = [nearest_idx3(57.5,F);nearest_idx3(62.5,F)];
    % De-Nanning
    powerCorrSort.nonansTotals = [];
    % Create NaNs index
    inds = find(squeeze(~isnan(powCorrTFR{1,ind1}.powspctrm(1,1,1,:)))); % Set D as all non-NaN numbers
    for z = 1:size(powCorrTFR{1,ind1}.powspctrm,1);
        for c = 1:length(powCorrTFR{1,ind1}.label);
            for j = 1:length(powCorrTFR{1,ind1}.powspctrm(1,1,:,1));
                % Fill with valid entries only
                powerCorrSort.nonansTotals(z,c,j,:) = powCorrTFR{1,ind1}.powspctrm(z,c,j,inds);
                10*log10(powCorrTFR{1,ind1}.powspctrm(z,c,j,inds));
                %powerCorrSort.nonansTotals(:,:,notchInd(1):notchInd(2),:) = NaN;
            end
        end
    end
    
    for z = 1:size(powerCorrSort.nonansTotals,1);
        for c = 1:length(powCorrTFR{1,ind1}.label);
            for t = 1:length(powerCorrSort.nonansTotals(1,1,1,:));
                % Fill with valid entries only
                powerCorrSort.nonansTotals(z,c,notchInd(1):notchInd(2),t) = NaN;
                powerCorrSort.nonansTotals(z,c,notchInd(1):notchInd(2),t) = interp1(find(~isnan(powerCorrSort.nonansTotals(z,c,:,t))),sq(powerCorrSort.nonansTotals(z,c,squeeze(~isnan(powerCorrSort.nonansTotals(z,c,:,t))),t)),find(isnan(powerCorrSort.nonansTotals(z,c,:,t))),'linear');
            end
        end
    end
    
    powerCorrSort.nonansBands1 = powerCorrSort.nonansTotals(:,:,(bandInd1:bandInd2),:);
    powerCorrSort.nonansBands2 = powerCorrSort.nonansTotals(:,:,(bandInd3:bandInd4),:);
    % Sorting Program
    
    for z = 1:size(powCorrTFR{1,ind1}.powspctrm,1);
        for c = 1:length(powCorrTFR{1,ind1}.label);
            for j = 1:2;
                % Repeat loop to acquire power for each time window
                for t = 1:size(powerCorrSort.nonansTotals,4)
                    % Master total frequency sort code line
                    powerCorrSort.mastertotals(z,c,:,t) = trapz(powerCorrSort.nonansTotals(z,c,freqRange,t));%-trapz(powerCorrSort.nonansTotals(z,c,(notchInd(1):notchInd(2)),t));
                    % First frequency band sort code line
                    powerCorrSort.masterbands(z,c,1,t) = trapz(powerCorrSort.nonansBands1(z,c,:,t));
                    % Second frequency band sort code line
                    powerCorrSort.masterbands(z,c,2,t) = trapz(powerCorrSort.nonansBands2(z,c,:,t));
                    if sum(bandInd1:bandInd2==31)==1
                        powerCorrSort.masterbands(z,c,1,t) = trapz(powerCorrSort.nonansBands1(z,c,:,t));%-trapz(powerCorrSort.nonansTotals(z,c,(notchInd(1):notchInd(2)),t)); %Notchout line
                    end
                    if sum(bandInd3:bandInd4==31)==1
                        powerCorrSort.masterbands(z,c,2,t) = trapz(powerCorrSort.nonansBands2(z,c,:,t));%-trapz(powerCorrSort.nonansTotals(z,c,(notchInd(1):notchInd(2)),t)); %Notchout line
                    end
                    powerCorrSort.masterpercents(z,c,j,t) = powerCorrSort.masterbands(z,c,j,t)/powerCorrSort.mastertotals(z,c,:,t); %Sort percent of total power for each frequency data point
                end
            end
        end
    end
    
    % Compute correlation between each frequency x channel
    numCmb = size(powerCorrSort.masterpercents,2)*size(powerCorrSort.masterpercents,3); %Create channel x frequency combination total
    p = reshape(powerCorrSort.masterpercents,[size(powerCorrSort.masterpercents,1), numCmb, size(powerCorrSort.masterpercents,4)]);
    %Use numCmb to create the nchoosek for all existing combinations
    varCmb = nchoosek(1:numCmb,2); 
    
    for z = 1:size(powerCorrSort.masterpercents,1);
        for e = 1:length(varCmb);
            [R] = corrcoef(p(z,varCmb(e,1),:),p(z,varCmb(e,2),:));
            powerCorrSort.mastercorr(z,e) = R(1,2);
        end
    end
    
    switch cfg.trialwindows
        case 'yes'
            powerCorrSort.trialWinString=[];
            for w = 1:size(p,2);
                p2 = p(:,w,:);
                p3 = sq(p2);
                p4 = reshape(p3,[(size(p3, 1)*size(p3, 2)), 1]);
                powerCorrSort.trialWinString(:,w) = p4(:,1);
            end
            for e = 1:length(varCmb);
                [R] = corrcoef(powerCorrSort.trialWinString(:,varCmb(e,1)),powerCorrSort.trialWinString(:,varCmb(e,2)));
                powerCorrSort.trialWinCorr(1,e) = R(1,2);
            end
        case 'no'
            TWCorr =[];
    end
    if iteration == 1
        powerCorrSort.mastercorrMeanPlot = zeros(size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1,size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1);
        powerCorrSort.mastercorrSTDPlot = zeros(size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1,size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1);
        powerCorrSort.mastercorrTWPlot = zeros(size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1,size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1);
    end
    m = mean(powerCorrSort.mastercorr,1)';
    s = std(powerCorrSort.mastercorr,0,1)';
    switch cfg.trialwindows        
        case 'yes'
            TW = (powerCorrSort.trialWinCorr)';
        case 'no'
            TWCorr =[];
    end   
    t = numCmb-1;
    z = numCmb;
    ii = numCmb-1;
    q = 1;
    c = 1;
    d = 2;
    b = numCmb;    
    if iteration == 1
        for j = 1:numCmb;
            powerCorrSort.mastercorrMeanPlot(d:b,j) = m(q:t,1);
            powerCorrSort.mastercorrStdPlot(d:b,j) = s(q:t,1);
            switch cfg.trialwindows
                case 'yes'
                    powerCorrSort.mastercorrTWPlot(d:b,j) = TW(q:t,1);
                case 'no'
                    TWCorr =[];
            end
            d = d+1;
            z = z-1;
            ii = ii-1;
            if j == 1;
                q = numCmb;
                t = numCmb+numCmb-3;
            end
            if j >= 2;
                q = q+z;
                t = (ii+t);
            end
        end
    end
    
    if iteration>=2
        if iteration==2
            d = 4*iteration-2;
            b = 4+iteration*4;
            c = (1+iteration*4)-4;
            g = c+6;
        end
        if iteration>2
            d = 4*iteration-2;
            b = 4+iteration*4;
            c = (1+iteration*4)-4;
            g = c+6;
        end
        for j = c:g;
            powerCorrSort.mastercorrMeanPlot(d:b,j) = m(q:t,1);
            powerCorrSort.mastercorrStdPlot(d:b,j) = s(q:t,1);
            switch cfg.trialwindows
                case 'yes'
                    powerCorrSort.mastercorrTWPlot(d:b,j) = TW(q:t,1);
                case 'no'
            end
            d = d+1;
            z = z-1;
            ii = ii-1;
            if j == c;
                q = numCmb;
                t = numCmb+numCmb-3;
            end
            if j > c;
                q = q+z;
                t = (ii+t);
            end
        end
    end
    % Plot Means
    if iteration == length(bands)-1;
        
        powerCorrSort.mastercorrMeanPlot(powerCorrSort.mastercorrMeanPlot==0) = NaN;
        %line(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'})
        
        MeanCorr = figure('units','normalized','position',[.7 .7 .8 .8]);
        pcolor(powerCorrSort.mastercorrMeanPlot)
        tall_str = sprintf(['\\fontsize{14}' blanks(1) '\\fontsize{11}']);
        h = title({'Average Channel x Frequency Correlations Across Trials' ;...
            [tall_str ]},...
            'FontWeight','bold',...
            'FontSize',11,...
            'FontName','Palatino Linotype');
        hc = colorbar('location','eastoutside','position',[.928 0.11 0.05 .815]);
        axis xy
        set(gca,'XTick',1.5:1:20,'XTickLabel',{'1','2','3','4'})
        set(gca,'YTick',1.5:1:20,'YTickLabel',{'1','2','3','4'})
        xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        % Set colormap
        colormap jet
        ax1 = gca;
        % Get position of first axes
        ax1_pos = ax1.Position; 
        ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');
        switch cfg.trialwindows
            case 'yes'
                powerCorrSort.mastercorrTWPlot(powerCorrSort.mastercorrTWPlot==0) = NaN;
                
                TWCorr = figure('units','normalized','position',[.7 .7 .8 .8]);
                pcolor(powerCorrSort.mastercorrTWPlot)
                tall_str = sprintf(['\\fontsize{14}' blanks(1) '\\fontsize{11}']);
                h = title({'Trial Window x Channel Frequency Correlations Across Trials' ;...
                    [tall_str ]},...
                    'FontWeight','bold',...
                    'FontSize',11,...
                    'FontName','Palatino Linotype');
                hc = colorbar('location','eastoutside','position',[.928 0.11 0.05 .815]);
                axis xy
                set(gca,'XTick',1.5:1:20,'XTickLabel',{'1','2','3','4'})
                set(gca,'YTick',1.5:1:20,'YTickLabel',{'1','2','3','4'})
                xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
                ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
                % Set colormap
                colormap jet
                ax1 = gca;
                % Get position of first axes
                ax1_pos = ax1.Position;
                ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');
            case 'no';
                TWCorr =[];
        end
        % Plot Real STDs
        powerCorrSort.mastercorrStdPlot(powerCorrSort.mastercorrStdPlot==0)=NaN;
                
        STDCorr = figure('units','normalized','position',[.7 .7 .8 .8]);
        pcolor(powerCorrSort.mastercorrStdPlot)
        tall_str = sprintf(['\\fontsize{14}' blanks(1) '\\fontsize{11}']);
        h = title({'Standard Deviations of Channel x Frequency Correlations Across Trials' ;...
            [tall_str ]},...
            'FontWeight','bold',...
            'FontSize',11,...
            'FontName','Palatino Linotype');
        hc = colorbar('location','eastoutside','position',[.928 0.11 0.05 .815]);
        axis xy
        set(gca,'XTick',1.5:1:20,'XTickLabel',{'1','2','3','4'})
        set(gca,'YTick',1.5:1:20,'YTickLabel',{'1','2','3','4'})
        xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        
        ax1=gca;
        % Get position of first axes
        ax1_pos = ax1.Position;
        ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');
    else
        MeanCorr = [];
        STDCorr = [];
        TWCorr = [];
    end
end
