it = [100:100:1000];
p = [100:50:1000];
for ii = 1:size(p,2)
    for jj = 1:size(it,2)
        for k = 1:it(jj)
            test(:,k) = randperm(22738,p(ii));
        end
        perc(ii,jj) = size(unique(reshape(test,1,size(test,1)*size(test,2))),2)/22738;
        clear test
    end

end
