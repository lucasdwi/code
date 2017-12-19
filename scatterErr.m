function h = scatterErr(x,y,err,fig,varargin)
%% Plot circles at (x,y) with error bars
%__________________________________________________________________________
% INPUTS
% x = x-axis values; format: integers, should be equal length vector as y
% y = y vlaues; format: vector of same length as x
% err = error values for error bars; format: vector of same length as x
% fig = whether to create new figure or not; format: either 1 or 0
% varargin
%   col = color of circles and lines
%   line = linestyle
%__________________________________________________________________________
% LLD 2017
%__________________________________________________________________________
%% Parse inputs
p = inputParser;
p.CaseSensitive = false;
addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric);
addRequired(p,'err',@isnumeric);
addRequired(p,'fig',@isnumeric);
addParameter(p,'col','k');
addParameter(p,'line','-');
addParameter(p,'mark','o');
addParameter(p,'markersize',36);
parse(p,x,y,err,fig,varargin{:});
%%
if fig == 1
    figure
end
if strcmpi(p.Results.mark,'line')
     % Get 1/4 of the difference between steps along x
    d = mean(diff(x)./4);
    for ii = 1:size(x,2)
        plot([x(ii)-d x(ii)+d],[y(ii) y(ii)],'LineStyle',p.Results.line,...
            'Color',p.Results.col)
    end
else
   h = scatter(x,y,p.Results.markersize,p.Results.col,'Marker',...
        p.Results.mark); 
end
hold on
for ii = 1:length(x)
   plot([x(ii) x(ii)],[y(ii)+err(ii) y(ii)-err(ii)],'LineStyle',...
       p.Results.line,'Color',p.Results.col) 
end