files = fileSearch('F:/irdmContinuous/','.mat');
for ii = 1:numel(files)
   s = dir(files{ii});
   procSize(ii) = s.bytes;
   if exist(['G:/GreenLab/data/irdmNew/split/',files{ii}(1:end-8),'.mat'],'file') == 2
       s = dir(['G:/GreenLab/data/irdmNew/split/',files{ii}(1:end-8),'.mat']);
       rawSize(ii) = s.bytes;
   else
       rawSize(ii) = NaN;
   end
end


figure
plot(rawSize,procSize,'o')

