for ci = 1:size(x,2)
   for ri = 1:size(x,1)
      z(ri,ci) = (x(ri,ci)-min(x(:,ci)))/(max(x(:,ci))-min(x(:,ci)));
   end
end
%%
figure
hold on
% LRN - Shell
plot(1,z(group(:,1)==1,1),'.k')
plot(1,z(group(:,1)==0,1),'.','Color',[0.5 0.5 0.5])
% LRN - Core
plot(1.5,z(group(:,3)==1,1),'.k')
plot(1.5,z(group(:,3)==0,1),'.','Color',[0.5 0.5 0.5])
% CPP - Shell
plot(2,z(group(:,1)==1,2),'.k')
plot(2,z(group(:,1)==0,2),'.','Color',[0.5 0.5 0.5])
% CPP - Core
plot(2.5,z(group(:,3)==1,2),'.k')
plot(2.5,z(group(:,3)==0,2),'.','Color',[0.5 0.5 0.5])
% Binge - Shell
plot(3,z(group(:,1)==1,3),'.k')
plot(3,z(group(:,1)==0,3),'.','Color',[0.5 0.5 0.5])
% Binge - Core
plot(3.5,z(group(:,3)==1,3),'.k')
plot(3.5,z(group(:,3)==0,3),'.','Color',[0.5 0.5 0.5])
% Rest - Shell
plot(4,z(group(:,1)==1,4),'.k')
plot(4,z(group(:,1)==0,4),'.','Color',[0.5 0.5 0.5])
% Rest - Core
plot(4.5,z(group(:,3)==1,4),'.k')
plot(4.5,z(group(:,3)==0,4),'.','Color',[0.5 0.5 0.5])
xlim([0.5 5])
title('Behavior Variation Across Shell and Core Responders')