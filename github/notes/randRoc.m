n = 20000;
p = 6;
thresh = 0:.001:1;

randProb = randi(100,n,1)./100;
tar = zeros(n,1);
ind = randi(n,p,1);
tar(ind) = 1;

[tpr,fpr] = rocCurve(tar,randProb,thresh);