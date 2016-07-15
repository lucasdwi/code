function [eventInds,eventTs] = eventInd(eventTs)
%% Finds indices within eventTs that match behavior labels
% eventTs = eventTs variable from converted data file given by
%   ConvertPl2All_Files.m (Stott)
% Necessary because not all data files have the same indices of events in
% eventTs
indO = find(not(cellfun('isempty',strfind(eventTs.label,'Orientation'))));
indRs = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (Start)'))));
indRe = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (End)'))));
indAs = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (Start)'))));
indAe = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (End)'))));
indBs = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (Start)'))));
indBe = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (End)'))));
indRm = 11;
% Set up index matrix 
eventInds = [];
eventInds(1:8,1) = [indO,indRs,indRe,indAs,indAe,indBs,indBe,indRm];
% Create Rest Middle
eventTs.label{11} = 'Rest (Middle)';
eventTs.t{11} = eventTs.t{1,indRs} + (eventTs.t{1,indRe} - eventTs.t{1,indRs})/2;
% Create generic rest, binge, and approach markers
eventTs.label{12} = 'Binge';
eventTs.label{13} = 'Rest';
eventTs.label{14} = 'Approach';