files1 = fileSearch('G:\GreenLab\code','.m');
files2 = fileSearch('G:\GreenLab\code\discovery','.m');
files3 = fileSearch('G:\GreenLab\code\notes','.m');
files4 = fileSearch('G:\GreenLab\code\old','.m');
files5 = fileSearch('G:\GreenLab\code\sp_package','.m');

files = [files1,files3,files4,files5];
for fi = 1:numel(files)
fid = fopen(files{fi});
res={};
while ~feof(fid)
  thisline = fgetl(fid);
  if ~ischar(thisline); break; end
  res{end+1,1} = thisline;
end
fclose(fid);
number_of_lines(fi) = numel(res);
end