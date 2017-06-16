function scatterErr(x,y,err,fig,varargin)
if fig == 1
    figure
end
if isempty(varargin)
    col = 'k';
else
    col = varargin{1};
end
hold on
plot(x,y,'o','Color',col)
for ii = 1:length(x)
   plot([x(ii) x(ii)],[y(ii)+err(ii) y(ii)-err(ii)],'-','Color',col) 
end