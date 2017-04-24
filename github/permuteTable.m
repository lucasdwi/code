function [] = permuteTable(pow1,pow2,coh1,coh2,powCorr1,powCorr2,ind,k)
% Get matrix of all combinations of rows of indices to keep, leaveing k out
inMat = nchoosek(1:size(ind,1),size(ind,1)-k)';
% Get matrix of all combinations of k rows of indices to get rid of and
% flip to match inMat
outMat = flipud(nchoosek(1:size(ind,1),k))';
% Preallocate
indIn = [];
indOut = [];

for ii = 1:size(inMat,2)
    indIn(:,ii) = sort(reshape(ind(inMat(:,ii),:),[numel(ind(inMat(:,ii),:)),1]));
    indOut(:,ii) = sort(reshape(ind(outMat(:,ii),:),[numel(ind(outMat(:,ii),:)),1]));
end
% Create tables of data corresponding to all possible combinations
for fsi = 1:size(pow1,1)
    T1{1,fsi} = [];
    for fi = 1:size(pow1,2)
        for ii = 1:size(indIn,2)
            chunk1 = [zeros(size(pow1{fsi,fi},3),1),reshape(pow1{fsi,fi}(:,:,:),[20,size(pow1{fsi,fi},3)])',reshape(coh1{fsi,fi}(:,:,:),[30,size(coh1{fsi,fi},3)])',powCorr1{fsi,fi}];
            chunk2 = [ones(size(pow2{fsi,fi},3),1),reshape(pow2{fsi,fi}(:,:,:),[20,size(pow2{fsi,fi},3)])',reshape(coh2{fsi,fi}(:,:,:),[30,size(coh2{fsi,fi},3)])',powCorr2{fsi,fi}];
            T1{ii,fsi} = vertcat(T1{ii,fsi},chunk1,chunk2);
        end
    end
end