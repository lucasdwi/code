function [cfg] = lassoNetCfg(naive,rand,normalize,foldGen,cvIterations,minTerm)
cfg.naive = naive;
cfg.rand = rand;
cfg.normalize = normalize;
cfg.foldGen = foldGen;
cfg.cvIterations = cvIterations;
cfg.minTerm = minTerm;
