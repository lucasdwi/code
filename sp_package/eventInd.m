function [eventInds] = eventInd(eventTs,eoi)
%% Finds indices within eventTs that match given behavior labels. This is 
% necessary because some data files may use a different order for event
% labels in eventTs - this depends on how the evenTs structure was created
% during video scoring. Works for both interval behaviors and scalars. In
% the case of scalar, the index will be repeated as if a start and stop
% exist in order to maintain matrix dimensions.
%
% N.B. given behavior labels (1) do NOT need to use the same case and (2)
% can be the start of a label (i.e. 'rest' will find 'resting'; will only
% compare the first N characters of labels in eventTs.label where N = the
% number of characters given in 'eoi' structure
%__________________________________________________________________________
% INPUTS:
% eventTs = eventTs variable from converted data file given by
%   ConvertPl2All_Files.m
% eoi = events of interest to be found; format = string cell array
%__________________________________________________________________________
% OUTPUTS: 
% eventInds = indices corresponding to eois; in the same order as eoi
%__________________________________________________________________________
% USE: 
% [eventInds] = eventInd(eventTs,{'binge',[0 5];'rest',[0 5]}) 
% Will extract the proper indices for the start and stop times of 'binge'
% and 'rest'; will work even if the rest behavior happens to be labelled
% 'resting'.
%__________________________________________________________________________
% DEPENDENCIES:
% logicFind.m
%__________________________________________________________________________
% LLD 2016-17
%% Go through eoi labels
eventInds = zeros(size(eoi,1),2);
for ei = 1:size(eoi,1)
   % Search for eoi label in eventTs.label
   inds = logicFind(1,strncmpi(eventTs.label,eoi{ei,1},...
       length(eoi{ei,1})),'==');
   if isempty(inds)
      error(['Warning: Behavior ',eoi{ei,1},' could not be found in '...
          'eventTs.label; make sure that spelling is correct in both '... 
          'eoi and eventTs.label.']) 
   end
   % If 2, assume that one is start and the other is end
   if length(inds) >2
       error([eoi{ei,1},' has a similarly named event;'...
           ' this will probably fuck up everything.'])
   elseif length(inds) == 2
      % Check if first label has 'Start' in it, if not make sure it has
      % 'End'
       if ~isempty(cell2mat(strfind(eventTs.label(inds(1)),'Start')))
           % Add indices to eventInds in correct order
           eventInds(ei,1:2) = inds(1:2);
       elseif ~isempty(cell2mat(strfind(eventTs.label(inds(1)),'End')))
               % Add indices to eventInds in correct order
               eventInds(ei,1:2) = inds(1:2)';
       end
   % Otherwise, only one occurence, concatenate to eventInds
   else
       eventInds(ei,1:2) = [inds(1),[]];
   end
end