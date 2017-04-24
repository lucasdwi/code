function [out] = shuffleMat(data,dim)
%% Shuffles matrix along dimension given
% If no dimension given, then shuffles all data
out = data(randperm(numel(data)));
for 
thisPerm = randperm(size(data,dim));
