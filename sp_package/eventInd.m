function [eventInds,eventTs,eventLabel,markers] = eventInd(eventTs,eoi)
%% Finds indices within eventTs that match behavior labels
% eventTs = eventTs variable from converted data file given by
%   ConvertPl2All_Files.m (Stott) 
% Necessary because not all data files have the same indices of events in
% eventTs

% Search label structure for indices
indO = find(not(cellfun('isempty',strfind(eventTs.label,'Orientation'))));
indRs = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (Start)'))));
indRe = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (End)'))));
indAs = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (Start)'))));
indAe = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (End)'))));
indBs = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (Start)'))));
indBe = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (End)'))));
indRm = 11;
%% Set up new index matrix 
eventInds = [];
eventInds(1:8,1) = [indO,indRs,indRe,indAs,indAe,indBs,indBe,indRm];
% Create 'Rest Middle' label
eventTs.label{11} = 'Rest (Middle)';
eventTs.t{11} = eventTs.t{1,indRs} + (eventTs.t{1,indRe} - eventTs.t{1,indRs})/2;
% Create generic rest, binge, and approach labels
eventTs.label{12} = 'Binge';
eventTs.label{13} = 'Rest';
eventTs.label{14} = 'Approach';
% Set up 'markers' vector
markers = [];
behaviors = {'app';'binge';'rest'};
behavInd = [eventInds(4), eventInds(6), eventInds(2)]; %[Approach, Binge, Rest]; gives eventTs index of behavior
markers(1,1) = find(not(cellfun('isempty',strfind(behaviors,eoi(1)))));
markers(1,2) = behavInd(markers(1,1));
if size(eoi,1) == 2
    markers(2,1) = find(not(cellfun('isempty',strfind(behaviors,eoi(2)))));
    markers(2,2) = behavInd(markers(2,1));
end
eventLabel = {eventTs.label{14},eventTs.label{12},eventTs.label{13}};