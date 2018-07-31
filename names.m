function [nameVect] = names(chan,band)
%% Generates a cell array of feature names in same order as CollateData.mat
%__________________________________________________________________________
% INPUTS
% chan = channel names; format = cell array
% band = frequency bands; format = cell array
%__________________________________________________________________________
% OUTPUTS
% nameVect = cell vector of all feature names
%__________________________________________________________________________
% LLD 2017
%% Power Names
for ii = 1:length(chan)
    for jj = 1:length(band)
        thisName = strcat(chan{ii},band{jj});
        cbNames{jj+((ii-1)*length(band))} = thisName;
    end
end
%% Coherence Names
cmbs = nchoosek(1:length(chan),2);
for iC = 1:size(cmbs,1)
    chlComb{iC} = strcat(chan(cmbs(iC,1)),chan(cmbs(iC,2)));
end
for ii = 1:length(chlComb)
    for jj = 1:length(band)
        thisName = strcat(chlComb{ii},band{jj});
        pbNames(jj+((ii-1)*length(band))) = thisName;
    end
end
%% Power Correlation Names
cmbs = nchoosek(1:length(chan),2);
c = 1;
for ii = 1:length(band)
   for jj = 1:size(cmbs,1)
        corrNames{c} = [chan{cmbs(jj,1)},band{ii},chan{cmbs(jj,2)}];
        c = c+1;
   end
end
%% Concatenate names
nameVect = [cbNames,pbNames,corrNames];