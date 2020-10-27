function med = circ_median(alpha,dim)
%
% med = circ_median(alpha)
%   Computes the median direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians, from 0 to 2*pi
%     [dim  compute along this dimension, default is 1, must 
%           be either 1 or 2 for circ_median]
%
%   Output:
%     mu		mean direction
%
%   circ_median can be slow for large datasets
%
% Update 2012
% PHB 3/19/2009
%
% References:
%   Biostatistical Analysis, J. H. Zar (26.6)
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 2
  dim = 1;
end

M = size(alpha);
med = NaN(M(3-dim),1);
for i=1:M(3-dim)
  if dim == 2
    beta = alpha(i,:)';
  elseif dim ==1
    beta = alpha(:,i);
  else
    error('circ_median only works along first two dimensions')
  end
  
  beta = mod(beta,2*pi);  %%Restrict range of angles between 0 and 2*pi
  
  %beta=sort(beta);
  n = size(beta,1);

  dd = circ_dist2(beta,beta);  % range between -pi and pi
  m1 = sum(dd>=0,1);  %%Angles greater than this angle
  m2 = sum(dd<=0,1);

  dm = abs(m1-m2);  %%%Minimum vallue means m1==m2, or median
  
  m = min(dm);
  idx = find(dm==m);
    

  if length(idx) > 1  %%Multiple points satisfy median criteria
      %disp('Ties detected.') %#ok<WNTAG>
      
      %%%Keep only points that lie within 90 degrees of each other
      %%%% (get rid of angles that are 180 degrees away; deal with them
      %%%% later)...
      delta_beta=abs(circ_dist2(beta(idx),beta(idx(1))));
      idx=idx(delta_beta<pi/2);
      % if length(delta_beta)>2
      %     fprintf('Many possible medians with delta: %s\n',mat2str(delta_beta',3));
      % end
      
  end


  md = circ_mean(beta(idx));  %Tentative median

  %%If median is 180 degrees off from center of mass, add 180 degrees
  circc_mean=circ_mean(beta);
  if abs(circ_dist(circc_mean,md)) > abs(circ_dist(circc_mean,md+pi))
    md = mod(md+pi,2*pi);  %Restrict to 0 to 2*pi
  end
  
  med(i) = md;
end

if dim == 2
  med = med';
end