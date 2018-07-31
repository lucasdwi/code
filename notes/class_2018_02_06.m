load('AH1_17_07_19_scored.mat')
[sdir,file,filter,dsf,thresh,onset,offset,foi,bands,overlap,cohMethod,eoi,vis,saveParent] = scbParamsMulti('AH1_17_07_19_scored.mat');
figure
plot(LFPTs.tvec,LFPTs.data(1,:))
[LFPTs.data] = filter60(LFPTs,adfreq,0);
[LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq);
minInt = max(diff(cell2mat(eoi(:,2)),1,2));
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,offset,...
    minInt,adfreq);
%%
figure
plot(LFPTs.tvec,LFPTs.data(1,:))
%%
hold on
for ii = 1:size(eventTs.t{6},1)
   plot([eventTs.t{6}(ii) eventTs.t{7}(ii)],[1 1],'-k')
end
%%
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
for ii = 1:size(trls{1,1}.sampleinfo,1)
   plot([LFPTs.tvec(trls{1,1}.sampleinfo(ii,1)) LFPTs.tvec(trls{1,1}.sampleinfo(ii,2))],[0.9 0.9]) 
end
%%
for ii = 1:size(trls{1,2}.sampleinfo,1)
   plot([LFPTs.tvec(trls{1,2}.sampleinfo(ii,1)) LFPTs.tvec(trls{1,2}.sampleinfo(ii,2))],[1.1 1.1],'-g') 
end
%%
inds = [nearest_idx3(1989,LFPTs.tvec):nearest_idx3(1990,LFPTs.tvec)];
figure
plot(LFPTs.tvec(inds),LFPTs.data(1,inds))
x = sigGen(adfreq,1,8,0.2,50);
hold on
plot(LFPTs.tvec(inds),x)
%%
inds2 = [nearest_idx3(1989.45,LFPTs.tvec):nearest_idx3(1989.6,LFPTs.tvec)];
y = sigGen(adfreq,0.15,77,0.07,220);
y = y +0.15;
plot(LFPTs.tvec(inds2),y)
xlim([1989.45 1989.6])
ylim([0.075 0.245])
