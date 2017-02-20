function [eventInds,eventTs,eventLabel] = eventInd(eventTs,eoi)
%% Finds indices within eventTs that match behavior labels
% Necessary because not all data files have the same indices of events in
% eventTs
% INPUTS
% eventTs = eventTs variable from converted data file given by
%   ConvertPl2All_Files.m (Stott) 

% OUTPUTS 
% markers = matrix of event indices; format row = behavior; in the
%           case of interval behaviors 1st column is start index, 2nd is
%           end index; in the case of discrete behavior, first column is
%           index and second is empty
%% Go through eoi labels
eventInds = [];
for ei = 1:size(eoi,1)
   % Search for eoi label in eventTs.label
   inds = logicFind(1,strncmpi(eventTs.label,eoi{ei,1},length(eoi{ei,1})),'==');
   % If 2, assume that one is start and other end
   if length(inds) == 2
      % Find which one corresponds to 
      % Check if first label has 'Start' in it, if not make sure it has
      % 'End'
       if ~isempty(cell2mat(strfind(eventTs.label(inds(1)),'Start')))
           % Add indices to eventInds in correct order
           eventInds(ei,1:2) = inds(1:2);
       else if ~isempty(cell2mat(strfind(eventTs.label(inds(1)),'End')))
               % Add indices to eventInds in correct order
               eventInds(ei,1:2) = inds(1:2)';
%                eventInds = [eventInds,inds(2),inds(1)];
           end
       end
   % Otherwise, only one occurence, concatenate to eventInds
   else
       eventInds(ei,1:2) = [inds(1),[]];
   end
end
% Autogenerate event labels from eoi
eventLabel = eoi(:,1)';