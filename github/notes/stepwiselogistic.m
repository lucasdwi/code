function stepwiselogistic(x,y)
[b,dev,stats] = glmfit(x,y(:,1),'binomial','link','logit');

maxdev = chi2inv(.95,1);
opt = statset('display','iter','TolFun',maxdev,'TolTypeFun','abs');
[in,history] = sequentialfs(@fitter2,x,y,'cv','none','nullmodel',true,'opt',opt,'direction','f');
dev = history.Crit;
nfeatures = sum(in);
plot((0:nfeatures)',dev,'b-x', nfeatures,dev(nfeatures+1),'ro')
[b,dev,st] = glmfit(x(:,in),y,'binomial');
% B = zeros(2,m+1);
B(1,[true,in]) = b';
B(2,[true,in]) = st.se';

function dev = fitter2(x,y)
[b,dev] = glmfit(x,y,'binomial');
