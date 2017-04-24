%% Examples of turning several indices into one continuous index vector
%% Two index counter
I = [1:n1];
J = [1:n2];
q = {};
for i = 1:length(I)
    for j = 1:length(J)
        qi = j+((i-1)*length(J));
        q = horzcat(q,qi);
    end
end
%% Three index counter
I = [1:n1];
J = [1:n2];
K = [1:n3];
q = [];
for i = 1:length(I)
    for j = 1:length(J)
        for k = 1:length(K)
            qi = (k+(j-1)*length(K))+(length(J)*length(K)*(i-1));
            q = horzcat(q,qi);
        end
    end
end