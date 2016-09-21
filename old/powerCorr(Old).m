function [STDCorr,MeanCorr,powerCorr]=powerCorr(TFRs,foi)
%% Use on ft_freqanalysis output to sort frequencies per trial per channel
%% Initialize color map
set(groot,'DefaultFigureColormap',jet);
%% Initialize frequency ranges and notch range
tic
powerCorr=[];
F=[];

bands.thet.limit = [4,7]; bands.alph.limit = [8,13]; bands.bet.limit = [15 30]; bands.lgam.limit = [45,65]; bands.hgam.limit = [70,90];
bandNames = fieldnames(bands);
F=TFRs{1, 1}.freq;
% n=1;
% c=1;
% z=1;
% j=1;
% t=1;
% ii=1;
% io=1;
%w=[2:45]; Hardcoded freqRange

freqRange = [bands.thet.limit(1)/foi(2):bands.hgam.limit(2)/foi(2)];
notchInd = [nearest_idx3(59.5,F);nearest_idx3(61.5,F)];
toc
%% Sorting Program
tic
for z = 1:length(TFRs{1,1}.powspctrm);    
    for c = 1:length(TFRs{1,1}.label);
        for j = 1:length(bandNames);
        % Create frequency band intervals from indices in F
        bands.(bandNames{j}).ind = [find(F>=bands.(bandNames{j}).limit(1),1),find(F<=bands.(bandNames{j}).limit(2),1,'last')];          
            for t = 1:(length(TFRs{1,1}.time))-101; %Repeat loop to acquire power for each time window
                for ii = length(bands.(bandNames{j}).ind(1):bands.(bandNames{j}).ind(2)), io=length(TFRs{1,1}.powspctrm(z,c,1:length(freqRange),(50+t)));
                    id = [bands.(bandNames{j}).ind(1):bands.(bandNames{j}).ind(2)]; %ni=[2:notchInd(1) notchInd(2):45]; jj=[bands.(bandNames{j}).ind(1):notchInd(1) notchInd(2):bands.(bandNames{j}).ind(2)];                         
                    [powerCorr.mastertotals(z,c,:,t)] = trapz(TFRs{1,1}.powspctrm(z,c,freqRange,(50+t)))-trapz(TFRs{1,1}.powspctrm(z,c,(notchInd(1):notchInd(2)),(50+t)));
                    [powerCorr.masterbands(z,c,j,t)] = trapz(TFRs{1,1}.powspctrm(z,c,id,(50+t))); %Master sort code line
                    if j == 4                            
                        [powerCorr.masterbands(z,c,j,t)] = (trapz(TFRs{1,1}.powspctrm(z,c,id,(50+t)))-trapz(TFRs{1,1}.powspctrm(z,c,(notchInd(1):notchInd(2)),(50+t))));
                    end
                end 
                [powerCorr.masterpercents(z,c,j,t)] = powerCorr.masterbands(z,c,j,t)/powerCorr.mastertotals(z,c,:,t);
            end
        end
    end
end 
toc
%% Compute correlation between each frequency x channel

powerCorr.mastercorr=[];
numCmb = size(powerCorr.masterpercents,2)*size(powerCorr.masterpercents,3);
p = reshape(powerCorr.masterpercents,[size(powerCorr.masterpercents,1), numCmb, size(powerCorr.masterpercents,4)]);
varCmb = nchoosek(1:numCmb,2);
tic
for z = 1:length(TFRs{1,1}.powspctrm);    
    for e = 1:length(varCmb); 
        [R] = corrcoef(p(z,varCmb(e,1),:),p(z,varCmb(e,2),:));
        powerCorr.mastercorr(z,e) = R(1,2);
    end    
end
toc
%% Reshape Correlogram Means
tic
powerCorr.mastercorrMeanPlot = [];
m = mean(powerCorr.mastercorr,1)';
s = std(powerCorr.mastercorr,0,1)';

d = 2;
t = numCmb-1;
z = numCmb;
ii = numCmb-1;
q = 1; 

for j = 1:numCmb;
    powerCorr.mastercorrMeanPlot(d:numCmb,j) = m(q:t,1);
    powerCorr.mastercorrStdPlot(d:numCmb,j) = s(q:t,1);
    d = d+1;
    z = z-1;
    ii = ii-1;
    if j == 1;
        q = 20;
        t = 37;
    end
    if j >= 2;
        q = q+z;
        t = (ii+t);
    end
    powerCorr.mastercorrMeanPlot(d-2,j) = 0;
    powerCorr.mastercorrStdPlot(d-2,j) = 0;
end    

powerCorr.mastercorrMeanPlot(numCmb+1,1:numCmb+1) = 0;
powerCorr.mastercorrStdPlot(numCmb+1,1:numCmb+1) = 0;
toc
%% Plot Means  

powerCorr.mastercorrMeanPlot(powerCorr.mastercorrMeanPlot==0) = NaN;
%line(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'})

MeanCorr = figure('units','normalized','position',[.7 .7 .8 .8]);
pcolor(powerCorr.mastercorrMeanPlot)
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
set(gca,'XTick',1.5:1:20,'XTickLabel',{'1','2','3','4'})
set(gca,'YTick',1.5:1:20,'YTickLabel',{'1','2','3','4'})
xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
%set(gca,'XTick',1:5:10)
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');


%% Plot Correlogram
% 
% powerCorr.mastercorrSTDPlot = [];
% 
% X = powerCorr.mastercorr;
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
%     powerCorr.mastercorrSTDPlot(d:numCmb,j)=m(q:t,1);
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
%     powerCorr.mastercorrSTDPlot(d-2,j)=0;
%     
% 
%   end    
%   powerCorr.mastercorrMeanPlot(21,1:21)=0;
%   
%   toc

%% Plot Real STDs  

  powerCorr.mastercorrStdPlot(powerCorr.mastercorrStdPlot==0)=NaN;
  

STDCorr=figure('units','normalized','position',[.7 .7 .8 .8]);
pcolor(powerCorr.mastercorrStdPlot)
tall_str = sprintf(['\\fontsize{14}' blanks(1) '\\fontsize{11}']);
h = title({'Standard Deviations of Channel x Frequency Correlations Across Trials' ;...
        [tall_str ]},...
  'FontWeight','bold',...
  'FontSize',11,...
  'FontName','Palatino Linotype');
%title('correlation')
%colorbar [5 1 0.8 0.9];
hc=colorbar('location','eastoutside','position',[.928 0.11 0.05 .815]);
%set(hc,'yaxisloc','top');
axis xy 
set(gca,'XTick',1.5:1:20,'XTickLabel',{'1','2','3','4'})
set(gca,'YTick',1.5:1:20,'YTickLabel',{'1','2','3','4'})
xlabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');
ylabel('Channels','FontWeight','bold','FontSize',11,'FontName','Palatino Linotype');

ax1=gca;
ax1_pos = ax1.Position; % position of first axes
%set(gca,'XTick',1:5:10)
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','XTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'FontName','Palatino Linotype','YAxisLocation','right','YTickLabel',{'-','\theta','-','\alpha','-','\beta','-','l\gamma','-','h\gamma','-'},'Color','none');


