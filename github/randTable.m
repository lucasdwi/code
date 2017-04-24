n = 15; p = 5;
randT = [];
randT(:,1) = round(rand(n,1));
%randT(:,1) = rand(n,1);
randT(:,2:p+1) = randn(n,p);
T.Base = array2table(randT);