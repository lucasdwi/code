function [nameVect] = names(chan,band)
%% Power Names
for i = 1:length(chan)
    for j = 1:length(band)
        thisName = strcat(chan{i},'-',band{j});
        cbNames{j+((i-1)*length(band))} = thisName;
    end
end
%% Coherence Names
cmbs = nchoosek(1:length(chan),2);
for iC = 1:size(cmbs,1)
    test{iC} = strcat(chan(cmbs(iC,1)),chan(cmbs(iC,2)));
end
for i = 1:length(chlComb)
    for j = 1:length(band)
        thisName = strcat(chlComb{i},'-',band{j});
        pbNames{j+((i-1)*length(band))} = thisName;
    end
end
%% Power Correlation Names
for ii = 1:size(corrNameParts,1)
   corrNames{ii} = [corrNameParts{ii,1},corrNameParts{ii,2},'-',corrNameParts{ii,3},corrNameParts{ii,4}] ;
end
[corrNames] = powCorrNames(band,chan);
%% Concatenate names
nameVect = [cbNames,pbNames,corrNames'];
%%
% chan = {'sl','sr','cl','cr'};
% chlComb = {'slsr'; 'slcl'; 'slcr'; 'srcl';'nscr';'clcr'};
% band = {'d','t','a','b','lg','hg'};
% corrNameParts = {'sl','sr','t','t';'sl','cl','t','t';'sl','cr','t','t';'sl','sl','t','a';'sl','sr','t','a';'sl','cl','t',...
%     'a';'sl','cr','t','a';'sr','cl','t','t';'sr','cr','t','t';'sr','sl','t','a';'sr','sr','t','a';'sr','cl','t','a';'sr','cr','t',...
%     'a';'cl','cr','t','t';'cl','sl','t','a';'cl','sr','t','a';'cl','cl','t','a';'cl','cr','t','a';'cr','sl','t','a';'cr','sr','t',...
%     'a';'cr','cl','t','a';'cr','cr','t','a';'sl','sr','a','a';'sl','cl','a','a';'sl','cr','a','a';'sl','sl','a','b';'sl','sr','a',...
%     'b';'sl','cl','a','b';'sl','cr','a','b';'sr','cl','a','a';'sr','cr','a','a';'sr','sl','a','b';'sr','sr','a','b';'sr','cl','a',...
%     'b';'sr','cr','a','b';'cl','cr','a','a';'cl','sl','a','b';'cl','sr','a','b';'cl','cl','a','b';'cl','cr','a','b';'cr','sl','a',...
%     'b';'cr','sr','a','b';'cr','cl','a','b';'cr','cr','a','b';'sl','sr','b','b';'sl','cl','b','b';'sl','cr','b','b';'sl','sl','b',...
%     'lg';'sl','sr','b','lg';'sl','cl','b','lg';'sl','cr','b','lg';'sr','cl','b','b';'sr','cr','b','b';'sr','sl','b','lg';'sr',...
%     'sr','b','lg';'sr','cl','b','lg';'sr','cr','b','lg';'cl','cr','b','b';'cl','sl','b','lg';'cl','sr','b','lg';'cl','cl','b',...
%     'lg';'cl','cr','b','lg';'cr','sl','b','lg';'cr','sr','b','lg';'cr','cl','b','lg';'cr','cr','b','lg';'sl','sr','lg','lg';...
%     'sl','cl','lg','lg';'sl','cr','lg','lg';'sl','sl','lg','hg';'sl','sr','lg','hg';'sl','cl','lg','hg';'sl','cr','lg',...
%     'hg';'sr','cl','lg','lg';'sr','cr','lg','lg';'sr','sl','lg','hg';'sr','sr','lg','hg';'sr','cl','lg','hg';'sr','cr',...
%     'lg','hg';'cl','cr','lg','lg';'cl','sl','lg','hg';'cl','sr','lg','hg';'cl','cl','lg','hg';'cl','cr','lg','hg';'cr',...
%     'sl','lg','hg';'cr','sr','lg','hg';'cr','cl','lg','hg';'cr','cr','lg','hg';'sl','sr','hg','hg';'sl','cl','hg','hg';...
%     'sl','cr','hg','hg';'sr','cl','hg','hg';'sr','cr','hg','hg';'cl','cr','hg','hg'};
