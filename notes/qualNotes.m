[x,~] = sigGen(100,30,2,1,0);
x = x(randperm(length(x)));
sig = x;
[x,~] = sigGen(100,1,80,5,0);
[x2,~] = sigGen(100,1,3,1,0);
x2 = x2(randperm(length(x2)));
stim1 = repmat([x,x2],1,10);
stim2 = repmat([x,x2(randperm(length(x2)))],1,10);
stim3 = repmat([x,x2(randperm(length(x2)))],1,10);
wash = sig(1:1001);
final = [sig,stim1,wash(randperm(length(wash))),stim2,wash(randperm(length(wash))),stim3,wash(randperm(length(wash)))];

final = final+rand(1,size(final,2));
tvec = linspace(0,120,length(final));
plot(tvec,final,'k','LineWidth',1)
set(gca,'XTick',0:10:120)
box off
%%
[x,~] = sigGen(100,5,2,1,0);
x = x.*rand(1,length(x));
x = x(randperm(length(x)));
sig = x;
x = [];
[x,~] = sigGen(100,5,20,5,0);
[x2,~] = sigGen(100,1,3,1,0);
x2 = x2.*rand(1,length(x2));
stim = [];
for ii = 1:19
    stim = [stim,x,x2(randperm(length(x2)))];
end
final = [sig,stim];
tvec = linspace(0,119,length(final));
figure
plot(tvec,final,'LineWidth',1)
xlim([0 20])