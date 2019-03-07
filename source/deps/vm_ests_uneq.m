function [mux,Rbar,kappa] = vm_ests_uneq(x,options,flag)
%VM_ESTS_UNEQ Maximum likelihood estimates of Von Mises parameters from data.
%  [MUX,RBAR,KAPPA] = VM_ESTS_UNEQ(X,OPTIONS,FLAG) estimates the mean vector (MUX),
%  the length of the mean vector (RBAR), and the concentration parameter (KAPPA)
%  of the Von Mises distribution.  X is assumed to be an n-by-2 matrix, where 
%  each row represents an x,y coordinate pair, which defines a bearing from
%  origin (0,0).  OPTIONS is an argument for FMINBND, for instance as set by
%  OPTIMSET.  FLAG is a boolean: 0 indicates that only MUX will be calculated;
%  1 indicates that RBAR and KAPPA will also be calculated
%
%  Calculates mean vector length, and thus KAPPA, taking account of individual
%  vector lengths.

n1 = size(x,1);
lx = sqrt(sum(x.^2,2));                % Get length of each vector
idx = find(lx>0);                        % Get rid of 0-length vectors
n = length(idx);

if n<n1
  x = x(idx,:);
end  
r = sum(x);
mux = r/sum(lx);                       % Mean x,y coordinate
if flag
  mu = 180/pi*atan2(mux(1),mux(2));      % Mean angle (use y/x for math convention)
  Rbar = norm(mux);                      % Length of mean vector
  s = 180/pi*sqrt(2*(1-Rbar));           % Angular standard deviation
  if Rbar<=(1/sqrt(n))
    kappa = 0;
  else
    kappa = fminbnd('diffkr',eps,5e4,options,Rbar,n); % Find kappa 
  end
end
