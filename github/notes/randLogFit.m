ri = 1; ci = 1;
testDat = [round(rand(100,1)),randn(100,1000)];
for r = 10:100;
    for c = 10:1000;
        %if c >= r*2
            [b,dev] = glmfit(testDat(1:r,2:c+1),testDat(1:r,1),'binomial');
            devMat(ri,ci) = dev;
        %else
        %    devMat(ri,ci) = NaN;
        %end
        ci = ci+1;
    end
    ci = 1;
    ri = ri+1;
end
%%
low = min(min(devMat));
high = max(max(devMat));
logAx = exp(linspace(log(low),log(high),64));

