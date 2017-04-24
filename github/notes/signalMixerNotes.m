fs = 1000;
t = 0:1/fs:4;
% Random frequencies
% f = randi([2 50],3,5);

% Set frequencies
f = [1,5,12,23,47];
a = rand(size(f,1),size(f,2));
p = 2*pi*(rand(size(f,1),size(f,2)));

% gcurve = normcdf(t,2.5,1);
% gcurve = normpdf(t,2.5,1);
% gcurve = rand(1,length(t));
gcurve = ones(1,length(t));
% noise = rand(1,length(t));
x = zeros(size(f,1),size(f,2),length(t));
sumx = zeros(size(f,2),length(t));
for fr = 1:size(f,1)
    for fi = 1:size(f,2)
        x(fr,fi,:) = a(fr,fi)*sin(f(fr,fi)*pi*t+p(fr,fi));
%     plot(t,x(fi,:))
    end
    sumx(fr,:) = sq(sum(x(fr,:,:),2))'.*gcurve;
    figure
    plot(t,sumx(fr,:))
end
%%

%%
fr = 100;
dt = 1/1000;
nBins = 10;
spikes = rand(5,nBins) < fr*dt;
figure
plot(sum(spikes,1))
%%
lenLog = 15000;
y = lognpdf(0:lenLog,log(2000),1);
figure
plot(0:lenLog,y)
len = 100000;
blank = zeros(15,len);
for ii = 1:15
   start = randi([0 len-lenLog]);
   blank(ii,start:start+lenLog) = y;
end
figure
plot(sum(blank,1));