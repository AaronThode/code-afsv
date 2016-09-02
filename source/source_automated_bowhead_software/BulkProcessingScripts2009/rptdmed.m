function [xm,ym,ints] = rptdmed(angs,src,r)
%RPTDMED Repeated median estimator of location from bearings.
%   [XM,YM,INTS] = RPTDMED(ANGS,src) assumes SRC is an m-by-2 matrix of 
%   Cartesian coordinates representing the source or origin of each bearing.
%   ANGS is an m-by-1 vector of bearings in radians (mathematical convention), 
%   ie, bearings correspond 1-1 with the coordinate pairs in SRC.  R is the
%   arbitrary length (same units as for SRC) of line segments for determining
%   whether bearings intersect (if intersection occurs at distance greater than
%   R, it is ignored - considered a non-intersection).  XM and YM are the x-
%   and y- coordinates of location, based on the repeated median estimator
%   described by Lenth (1981).  INTS is a logical vector with the same
%   dimension as ANGS (1 = the corresponding bearing has at least one
%   intersection with another bearing; 0 = no intersections).  The repeated
%   median depends on determining all points of intersection.  First, the
%   component-wise median is calculated across all of the intersection points
%   for each bearing.  The estimated location is the median of these medians.
%
%   Lenth,R.V. 1981.  On finding the source of a signal.
%     Technometrics 23:149-154.

m = length(angs);
m1 = m-1;
one = ones(1,m1);
idx = [1:m];
ic = zeros(m);
xc = ic;  yc = ic;
for i = 1:m1,
  j = (i+1):m;
  one = ones(1,m-i);
  x1 = src(i,1);           % Origins of bearings
  y1 = src(i,2);
  x2 = src(j,1)';
  y2 = src(j,2)';
  x1 = [x1; x1 + r*cos(angs(i,:))] * one;  % Endpoints of segments
  y1 = [y1; y1 + r*sin(angs(i,:))] * one;
  x2 = [x2; x2 + r*cos(angs(j,:))'];
  y2 = [y2; y2 + r*sin(angs(j,:))'];
  ic(i,j) = iscross(x1,y1,x2,y2);       % Determine whether segments intersect
  ic2 = find(ic(i,j));
  j2 = j(ic2);
  if ~isempty(ic2),                     % Intersection points
    [xc(i,j2),yc(i,j2)] = intsecl(x1(:,ic2),y1(:,ic2),x2(:,ic2),y2(:,ic2));
  end
end
% ic = ic + ic' + diag(NaN*ones(1,m));
ic = ic + ic';
ints = any(ic,2);
xc = xc + xc';
yc = yc + yc';
xc(~xc) = NaN;
yc(~yc) = NaN;
xm = nanmedd(nanmedd(xc,1),2);
ym = nanmedd(nanmedd(yc,1),2);