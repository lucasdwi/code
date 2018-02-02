function [r,rVect] = powerCorr(psdTrls)
%% Calculates power correlations within frequency bands across channels
% INPUTS:
% psdTrls = power-spectra for each behavior of interest; format: structure
%   created by powerComp.m where each behavior is in its own cell with
%   .relPow as a field containing the relative power in each band for every
%   trial
%__________________________________________________________________________
% OUTPUTS:
% r = cell array of correlation values; format: cell array with 3D arrays
%   containing all correlation values (channel,channel,band)
% rVect = vectorized version of r; because of intrinsic symmetry of
%   correlation matrices, takes just lower half of matrix and puts into row
%   vector. For example, if 4 channels then the vector of r values will be
%   for the following correlations: 1-2, 1-3, 1-4, 2-3, 2-4, and 3-4.
%__________________________________________________________________________
%% LLD and MAC 2017
%% Initialize
% Find empty events
empt = cellfun(@isempty,psdTrls);
%% Preallocate
% Get number of bands from .relPow field in first psdTrls cell that isn't
% empty
nBands = size(psdTrls{1,logicFind(0,empt,'==','first')}.relPow,1);
nEvents = size(psdTrls,2);
r = cell(1,nEvents);
p = cell(1,nEvents);
rVect = cell(1,nEvents);
%% Go through each event and all bands getting correlations of all channels
% to eachother
for iE = logicFind(0,empt,'==')
    if size(psdTrls{1,iE}.Pow,3) > 1
        for iB = 1:nBands
            % Calculate correlations
            [r{1,iE}(:,:,iB),p{1,iE}(:,:,iB)] = corrcoef(squeeze(psdTrls{1,iE}.relPow(iB,:,:))');
            % Extract just lower half of correlation matrix
            thisR = squeeze(r{1,iE}(:,:,iB));
            % Place in rVect which will be ordered by column then row (i.e. all
            % rows of column one, then 2, etc.)
            rVect{1,iE}(:,iB) = thisR(logicFind(0,triu(thisR),'=='));
        end
    else
        r{1,iE} = NaN;
        rVect{1,iE} = NaN;
    end
end
