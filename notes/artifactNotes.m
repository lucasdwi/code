load('ZA1_NAcCont250A_2018_04_13.mat')
%% Use eventTs.t{1} and LFPTs.data(5,:)
delay = 9;
artifact = zeros(size(eventTs.t{1},1),delay);

for ii = 1:size(eventTs.t{1},1)
    x = nearest_idx3(eventTs.t{1}(ii),LFPTs.tvec);
    y = nearest_idx3(eventTs.t{1}(ii),LFPTs.tvec)+delay-1;
    artifact(ii,:) = LFPTs.data(5,x:y);
end

mArtifact = mean(artifact,1);
sArtifact = std(artifact,[],1);
figure
shadedErrorBar(0:0.001:0.008,mArtifact,sArtifact)
q1 = mean(artifact(1:175584,:),1);
q2 = mean(artifact(175585:351168,:),1);
q3 = mean(artifact(351169:526752,:),1);
q4 = mean(artifact(526753:702336,:),1);
hold on
plot(0:0.001:0.008,q1)
plot(0:0.001:0.008,q2)
plot(0:0.001:0.008,q3)
plot(0:0.001:0.008,q4)
%%
test = LFPTs.data(5,:);
for ii = 1:size(eventTs.t{1},1)
    x = nearest_idx3(eventTs.t{1}(ii),LFPTs.tvec);
    y = nearest_idx3(eventTs.t{1}(ii),LFPTs.tvec)+delay-1;
    test(x:y) = test(x:y)-mArtifact;
end