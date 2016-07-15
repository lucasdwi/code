function fileCycle
%% Run spectcompbase.m
files = {'H10BaseSep27','H10FoodDepSep25','H13BaseSep27','H13FoodDepSep25','H14BaseSep29','H14FoodDepSep28','H15BaseSep29','H15FoodDepSep28','I1BaseNov9','I1FoodDepNov13','I2BaseNov9','I2FoodDep24Dec16','I3BaseNov9','I3FoodDep24Nov2','I4BaseSep24','I4FoodDepSep30','I6BaseSep18','I6FoodDepSep30','I8BaseSep23','I8FoodDepSep30','I11BaseOct30','I11FoodDep24Dec16','I12BaseNov12','I12FoodDepNov13'};
for i = 1:length(files)
    [LFPTsNaN,nNaN,indSkp,trl1,trl2,clnTrls,clnEvents,relPower,stdPower,psdTrls,TFR_event1,TFR_event2,fd1,fd2,avgCoh,relCoh,stdCoh] = spectcompbase(files{i},5,2.5,5,17000,3,1.5,2,3,[1 2 150],0.5,[0 0.005 3],[0 0.005 3]);
    close all; clearvars -except files i;
end
%% Extract all PSD and Coherence data
% if func == 2
%     % Set base and fooddep trials
%     
%     % Collect PSD and coherence data
%     [T,allVars] = tabulateData(files);
%     [models] = lineReg(T,allVars);
% end



    