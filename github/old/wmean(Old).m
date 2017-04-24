function [awx] = wmean(x,w)
% Uses vector of weights (w) corresponding to 3rd dimension of data (x)
% to compute a weighted average
for wi = 1:numel(w)
    for si = 1:size(x{wi},3)
        % Multiply each file PSD, coh, powCorr with weight matrix
        wx{wi}(:,:,si) = x{wi}(:,:,si)*w{wi}(si);
        % Sum across files
        swx{wi} = sum(wx{wi},3);
        % Divide by sum of weights
        awx{wi} = swx{wi}./sum(w{wi});
    end
end
