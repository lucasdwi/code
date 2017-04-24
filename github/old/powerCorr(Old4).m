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
    
    %iteration=4; %Tune for specific iterations if desired
    
    if iteration == 1      %Resets data struct for each new data set
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
    % z = 1;
    % c = 1;
    % j = 1;
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

    
    %
    
    
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
    
    varCmb = nchoosek(1:numCmb,2); %Use numCmb to create the nchoosek for all existing combinations
    
    for z = 1:size(powerCorrSort.masterpercents,1);
        for e = 1:length(varCmb);
            [R] = corrcoef(p(z,varCmb(e,1),:),p(z,varCmb(e,2),:));
            powerCorrSort.thisCorr(z,e) = R(1,2);
        end
    end
    
    
    if iteration == 1
        powerCorrSort.masterCorr(1:size(powerCorrSort.thisCorr,1),(1:size(powerCorrSort.thisCorr,2)-6)) = powerCorrSort.thisCorr(:,1:22);
    end
    
    if iteration >= 2
        powerCorrSort.masterCorr(1:size(powerCorrSort.thisCorr,1),((iteration-1)*size(powerCorrSort.thisCorr,2)-6*(iteration-1)+1):((iteration)*size(powerCorrSort.thisCorr,2)-6*(iteration))) = powerCorrSort.thisCorr(:,1:22);
    end
    
    if iteration == length(bands)-1
        powerCorrSort.masterCorr(1:size(powerCorrSort.thisCorr,1),(iteration-1)*size(powerCorrSort.thisCorr,2)-17:iteration*size(powerCorrSort.thisCorr,2)-18) = powerCorrSort.thisCorr(:,:);
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
            
            
            
            
            %for z = 1:size(powerCorrSort.masterpercents,1);
            for e = 1:length(varCmb);
                [R] = corrcoef(powerCorrSort.trialWinString(:,varCmb(e,1)),powerCorrSort.trialWinString(:,varCmb(e,2)));
                powerCorrSort.trialWinCorr(1,e) = R(1,2);
            end
            %end
            
            
        case 'no'
            
            TWCorr =[];
    end
    
    if iteration == 1
        powerCorrSort.masterCorrMeanPlot = zeros(size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1,size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1);
        powerCorrSort.masterCorrStdPlot = zeros(size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1,size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1);
        powerCorrSort.masterCorrTwPlot = zeros(size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1,size(bands,1)*size(powCorrTFR{1,ind1}.label,2)+1);
    end
    
    m = mean(powerCorrSort.thisCorr,1)';
    s = std(powerCorrSort.thisCorr,0,1)';
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
            powerCorrSort.masterCorrMeanPlot(d:b,j) = m(q:t,1);
            powerCorrSort.masterCorrStdPlot(d:b,j) = s(q:t,1);
            switch cfg.trialwindows
                case 'yes'
                    powerCorrSort.masterCorrTwPlot(d:b,j) = TW(q:t,1);
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
            %    powerCorrSort.masterCorrMeanPlot(d-2,j) = 0;
            %    powerCorrSort.masterCorrStdPlot(d-2,j) = 0;
        end
    end
    
    if iteration>=2
        if iteration>=2
            d = size(powCorrTFR{1,ind1}.label,2)*iteration-(size(powCorrTFR{1,ind1}.label,2)-2);
            b = size(powCorrTFR{1,ind1}.label,2)+iteration*size(powCorrTFR{1,ind1}.label,2);
            c = (1+iteration*size(powCorrTFR{1,ind1}.label,2))-size(powCorrTFR{1,ind1}.label,2);
            g = c+size(powCorrTFR{1,ind1}.label,2)+(size(powCorrTFR{1,ind1}.label,2)-2);
        end
        %    if iteration>2
        %        d = 4*iteration-2;
        %        b = 4+iteration*4;
        %        c = (1+iteration*4)-4;
        %        g = c+6;
        %    end
        for j = c:g;
            powerCorrSort.masterCorrMeanPlot(d:b,j) = m(q:t,1);
            powerCorrSort.masterCorrStdPlot(d:b,j) = s(q:t,1);
            switch cfg.trialwindows
                case 'yes'
                    powerCorrSort.masterCorrTwPlot(d:b,j) = TW(q:t,1);
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
    %powerCorrSort.masterCorrMeanPlot(numCmb+1,1:numCmb+1) = 0;
    %powerCorrSort.masterCorrStdPlot(numCmb+1,1:numCmb+1) = 0;
     
    
    
    % Plot Means
    
    
    if iteration == length(bands)-1;
        
        
        %set(groot,'DefaultFigureColormap',jet)
        
        powerCorrSort.masterCorrMeanPlot(powerCorrSort.masterCorrMeanPlot==0) = NaN;
        %line(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'})
        
        MeanCorr = figure('units','normalized','position',[.7 .7 .8 .8]);
        pcolor(powerCorrSort.masterCorrMeanPlot)
        tall_str = sprintf(['\\fontsize{14}' blanks(1) '\\fontsize{11}']);
        h = title({'Average Channel x Frequency Correlations Across Trials' ;...
            [tall_str ]},...
            'FontWeight','bold',...
            'FontSize',11,...
            'FontName','Palatino Linotype');
        %title('correlation')
        %colorbar [5 1 0.8 0.9];
        hc = colorbar('location','eastoutside','position',[.928 0.11 0.05 .815]);
        %set(hc,'yaxisloc','top');
        axis xy
        if size(powerCorrSort.masterCorrMeanPlot,1)>=9
            disp('Warning, too many channels, will not plot figure correctly.')
        end
        if size(powerCorrSort.masterCorrMeanPlot,1)<=2
            disp('Warning, too few channels, will not plot figure correctly.')
        end
        if  size(powCorrTFR{1,ind1}.label,2) == 3
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 4
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 5
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 6
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 7
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6','7'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6','7'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 8
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6','7','8'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6','7','8'})
            
        end
        
        
        xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        % Set colormap
        colormap jet
        ax1 = gca;
        ax1_pos = ax1.Position; % position of first axes
        %set(gca,'XTick',1:5:10)
        ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');
        
        
        switch cfg.trialwindows
            
            case 'yes'
                
                powerCorrSort.masterCorrTwPlot(powerCorrSort.masterCorrTwPlot==0) = NaN;
                
                TWCorr = figure('units','normalized','position',[.7 .7 .8 .8]);
                pcolor(powerCorrSort.masterCorrTwPlot)
                tall_str = sprintf(['\\fontsize{14}' blanks(1) '\\fontsize{11}']);
                h = title({'Trial Window x Channel Frequency Correlations Across Trials' ;...
                    [tall_str ]},...
                    'FontWeight','bold',...
                    'FontSize',11,...
                    'FontName','Palatino Linotype');
                %title('correlation')
                %colorbar [5 1 0.8 0.9];
                hc = colorbar('location','eastoutside','position',[.928 0.11 0.05 .815]);
                %set(hc,'yaxisloc','top');
                axis xy
                if  size(powCorrTFR{1,ind1}.label,2) == 3
                    set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3'})
                    set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3'})
                    
                elseif size(powCorrTFR{1,ind1}.label,2) == 4
                    set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4'})
                    set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4'})
                    
                elseif size(powCorrTFR{1,ind1}.label,2) == 5
                    set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5'})
                    set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5'})
                    
                elseif size(powCorrTFR{1,ind1}.label,2) == 6
                    set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6'})
                    set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6'})
                    
                elseif size(powCorrTFR{1,ind1}.label,2) == 7
                    set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6','7'})
                    set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6','7'})
                    
                elseif size(powCorrTFR{1,ind1}.label,2) == 8
                    set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6','7','8'})
                    set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6','7','8'})
                    
                end
                xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
                ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
                % Set colormap
                colormap jet
                ax1 = gca;
                ax1_pos = ax1.Position; % position of first axes
                %set(gca,'XTick',1:5:10)
                ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');
                
            case 'no';
                TWCorr =[];
        end
        
        
        %end
        
        
        % Plot Correlogram
        %
        % powerCorrSort.masterCorrSTDPlot = [];
        %
        % X = powerCorrSort.masterCorr;
        % m = std(X,0,1);
        % m = m';
        % N = 19;
        % numCmb=20;
        % q=1;
        % d=2;
        % j=1;
        % t=19;
        % z=20;
        % ii=19;
        %
        % %for C=1:length(M);
        %    %for D=1:length(20)
        %
        %    tic
        %   for j=1:20;
        %
        %     powerCorrSort.masterCorrSTDPlot(d:numCmb,j)=m(q:t,1);
        %
        %
        %     d=d+1;
        %
        %     z=z-1;
        %     ii=ii-1;
        %
        %     if j==1;
        %     q=20;
        %     t=37;
        %     end
        %
        %     if j>=2;
        %         q=q+z;
        %         t=(ii+t);
        %
        %
        %     end
        %
        %
        %
        %     powerCorrSort.masterCorrSTDPlot(d-2,j)=0;
        %
        %
        %   end
        %   powerCorrSort.masterCorrMeanPlot(21,1:21)=0;
        %
        %   toc
        
        % Plot Real STDs
        
        powerCorrSort.masterCorrStdPlot(powerCorrSort.masterCorrStdPlot==0)=NaN;
        
        
        STDCorr = figure('units','normalized','position',[.7 .7 .8 .8]);
        pcolor(powerCorrSort.masterCorrStdPlot)
        tall_str = sprintf(['\\fontsize{14}' blanks(1) '\\fontsize{11}']);
        h = title({'Standard Deviations of Channel x Frequency Correlations Across Trials' ;...
            [tall_str ]},...
            'FontWeight','bold',...
            'FontSize',11,...
            'FontName','Palatino Linotype');
        %title('correlation')
        %colorbar [5 1 0.8 0.9];
        hc = colorbar('location','eastoutside','position',[.928 0.11 0.05 .815]);
        %set(hc,'yaxisloc','top');
        axis xy
        if  size(powCorrTFR{1,ind1}.label,2) == 3
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 4
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 5
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 6
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 7
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6','7'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6','7'})
            
        elseif size(powCorrTFR{1,ind1}.label,2) == 8
            set(gca,'XTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'XTickLabel',{'1','2','3','4','5','6','7','8'})
            set(gca,'YTick',1.5:1:size(powerCorrSort.masterCorrMeanPlot,1),'YTickLabel',{'1','2','3','4','5','6','7','8'})
            
        end
        xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
        % Set colormap
        colormap jet
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        %set(gca,'XTick',1:5:10)
        ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');
        
        
        
        
    else
        
        MeanCorr = [];
        STDCorr = [];
        TWCorr = [];
        
        
        
    end
end
