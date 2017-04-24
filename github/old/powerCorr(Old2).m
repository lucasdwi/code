function [STDCorr,MeanCorr,powerCorrSort,freqRange,notchInd,numCmb] = powerCorr(TFRs,bands)
%% Use on ft_freqanalysis output to sort frequencies per trial per channel
%% Initialize color map
set(groot,'DefaultFigureColormap',jet);
%% Initialize frequency ranges and notch range
tic
powerCorrSort = [];

freqstart = find(TFRs{1,1}.freq>=bands{1,2}(1),1);
freqend = find(TFRs{1,1}.freq<=bands{end,2}(2),1,'last');
freqRange = freqstart:freqend; %Create indices for total frequencies between first bands to last bands, set at function level as 2Hz:90Hz

notchInd = [nearest_idx3(59.5,TFRs{1,1}.freq);nearest_idx3(61.5,TFRs{1,1}.freq)]; %Create notch index
toc
%% De-Nanning 

z=1;
c=1;
j=1;
N=1;
d=squeeze(isnan(TFRs{1,1}.powspctrm(z,c,j,:))); %Create NaNs index
D=find(N>d); %Set D as all real numbers

tic

for z = 1:length(TFRs{1,1}.powspctrm);   
    for c = 1:length(TFRs{1,1}.label);
        for j = 1:length(TFRs{1,1}.powspctrm(1,1,:,1));
            for t=1:(length(D));
                                                       
                                              
               powerCorrSort.nonans(z,c,j,:) = TFRs{1,1}.powspctrm(z,c,j,D); %Fill with valid entries only
                        
    
             end
        end
    end
end

toc

%% Sorting Program
tic
for z = 1:length(TFRs{1,1}.powspctrm);    
    for c = 1:length(TFRs{1,1}.label);
        for j = 1:length(bands);
            bandInd(j,1) = find(TFRs{1,1}.freq>=bands{j,2}(1),1);
            bandInd(j,2) = find(TFRs{1,1}.freq<=bands{j,2}(2),1,'last');
            
            for t = 1:(length(powerCorrSort.nonans(z,c,j,:))); %Repeat loop to acquire power for each time window
                 
                
       
                    powerCorrSort.mastertotals(z,c,:,t) = trapz(powerCorrSort.nonans(z,c,freqRange,(t)))-trapz(powerCorrSort.nonans(z,c,(notchInd(1):notchInd(2)),(t))); %Master total frequency sort code line
                    powerCorrSort.masterbands(z,c,j,t) = trapz(powerCorrSort.nonans(z,c,(bandInd(j,1):bandInd(j,2)),(t))); %Master frequency band sort code line
                    
                    if j == 4                            
                        powerCorrSort.masterbands(z,c,j,t) = (trapz(powerCorrSort.nonans(z,c,(bandInd(j,1):bandInd(j,2)),(t)))-trapz(powerCorrSort.nonans(z,c,(notchInd(1):notchInd(2)),(t)))); %Notchout line
                    end
        
                powerCorrSort.masterpercents(z,c,j,t) = powerCorrSort.masterbands(z,c,j,t)/powerCorrSort.mastertotals(z,c,:,t); %Sort percent of total power for each frequency data point
            end
        end
    end
end 
toc
%% Compute correlation between each frequency x channel


numCmb = size(powerCorrSort.masterpercents,2)*size(powerCorrSort.masterpercents,3); %Create channel x frequency combination total
p = reshape(powerCorrSort.masterpercents,[size(powerCorrSort.masterpercents,1), numCmb, size(powerCorrSort.masterpercents,4)]);
varCmb = nchoosek(1:numCmb,2); %Use numCmb to create the nchoosek for all existing combinations
tic
for z = 1:length(TFRs{1,1}.powspctrm);    
    for e = 1:length(varCmb); 
        [R] = corrcoef(p(z,varCmb(e,1),:),p(z,varCmb(e,2),:));
        powerCorrSort.mastercorr(z,e) = R(1,2);
    end    
end
toc
%% Reshape Correlogram Means
tic

m = mean(powerCorrSort.mastercorr,1)';
s = std(powerCorrSort.mastercorr,0,1)';

d = 2;
t = numCmb-1;
z = numCmb;
ii = numCmb-1;
q = 1; 

for j = 1:numCmb;
    powerCorrSort.mastercorrMeanPlot(d:numCmb,j) = m(q:t,1);
    powerCorrSort.mastercorrStdPlot(d:numCmb,j) = s(q:t,1);
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
    powerCorrSort.mastercorrMeanPlot(d-2,j) = 0;
    powerCorrSort.mastercorrStdPlot(d-2,j) = 0;
end    

powerCorrSort.mastercorrMeanPlot(numCmb+1,1:numCmb+1) = 0;
powerCorrSort.mastercorrStdPlot(numCmb+1,1:numCmb+1) = 0;
toc
%% Plot Means  

set(groot,'DefaultFigureColormap',jet)

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
% powerCorrSort.mastercorrSTDPlot = [];
% 
% X = powerCorrSort.mastercorr;
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
%     powerCorrSort.mastercorrSTDPlot(d:numCmb,j)=m(q:t,1);
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
%     powerCorrSort.mastercorrSTDPlot(d-2,j)=0;
%     
% 
%   end    
%   powerCorrSort.mastercorrMeanPlot(21,1:21)=0;
%   
%   toc

%% Plot Real STDs  

  powerCorrSort.mastercorrStdPlot(powerCorrSort.mastercorrStdPlot==0)=NaN;
  

STDCorr=figure('units','normalized','position',[.7 .7 .8 .8]);
pcolor(powerCorrSort.mastercorrStdPlot)
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


