%% Playing with coherence parameters: window size to accurately estimate 
% cohernce

fs = 500; dt = 1./fs;
t = [0 5]; tvec = t(1):dt:t(2)-dt;
 
f1 = 1:1:100; f2 = 1:1:100; phi = 0:0.01:2;
for fi1 = 1:100
    data1(fi1,:) = sin(2*pi*f1(fi1)*tvec)+0.1*randn(size(tvec));
end
for fi2 = 1:100
    for phiI = 1:201
        data2(fi2,:,phiI) = sin(2*pi*f2(fi2)*tvec+(pi*phi(phiI)))+0.1*randn(size(tvec)); % phase-shifted version of data1
    end
end
%%
c = 1;
for wi = 1:10
    for fi1 = 1:100
        [cxy(fi1,:,c),f] = mscohere(data1(fi1,:),data2(fi1,:,1),wi*100,ceil((wi*100)*0.75),1:100,fs);
    end
    c = c+1;
end