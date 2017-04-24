one = 15;
zero = 5;

test = [ones(1,one),zeros(1,zero)];
% Using 5-fold
cmbs = nchoosek(1:one+zero,length(test)/5);
both = 0;
for ii = 1:size(cmbs,1)
   uChk = unique(test(cmbs(ii,:))); 
   if size(uChk,2) == 2
       both = both + 1;
   end
end
perc = both/size(cmbs,1);
