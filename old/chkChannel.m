function chkChannel(LFPTsNaN,m)
%%
% Calculate root-mean squares for each channel
for c = 1:size(LFPTsNaN.data,1)
    r(c) = rms(LFPTsNaN.data(c,~isnan(LFPTsNaN.data(c,:))));
end
% Check if RMS is below specified value (min)
if any(r<m)
    disp('It looks like at least on of the channels might be bad; please review the plots, potential bad channel(s) in red. Press Ctl+C to exit or any other key to continue.')
    for c = 1:size(LFPTsNaN.data,1)
        figure;
        f{c} = plot(LFPTsNaN.data(c,:));
        if r(c) < m
            f{c}.Color = 'r';
        end       
    end
    pause
end
    