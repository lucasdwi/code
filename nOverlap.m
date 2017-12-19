function [xiny,yinx,percOver] = nOverlap(x,y)
%%
% Delaunay triangulation
xTri = delaunayn(x);
yTri = delaunayn(y);
% Looks for nearest enclosing simplex, if doesn't exit (i.e. point is
% outside) returns NaN
yinxTri = tsearchn(x,xTri,y);
xinyTri = tsearchn(y,yTri,x);
% Finds percent of points outside the other set
xiny = sum(~isnan(xinyTri))/size(x,1);
yinx = sum(~isnan(yinxTri))/size(y,1);
% Finds percent overlap of all points
percOver = (sum(~isnan(xinyTri))+sum(~isnan(yinxTri)))/(size(x,1)+size(y,1));