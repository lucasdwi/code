function [pick] = marchMadness(seed1,seed2)
% Subtract seed from 17 to account for lower numbers being better
seeds = [zeros(1,17-seed1),ones(1,17-seed2)];
% Add one s.t., if 0 seed1 and if 1 seed2
pick = seeds(randi(numel(seeds)))+1;

