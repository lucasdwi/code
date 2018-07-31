function [newX,newY] = randRoc(x,y)
sizes = cell2mat(cellfun(@size,x,'UniformOutput',0)');
samp = max(sizes(:,1));
for ii = 1:length(x)
   if isequal(x{ii},[0;1])
       newX(:,ii) = 0:1/(samp-1):1;
       newY(:,ii) = 0:1/(samp-1):1;
   else
       newX(:,ii) = x{ii};
       newY(:,ii) = y{ii};
   end
end