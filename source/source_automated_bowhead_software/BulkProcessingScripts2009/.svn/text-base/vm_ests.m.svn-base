function [mu,kappa,sd] = vm_ests(x,options)
%VM_ESTS Maximum likelihood estimates of Von Mises parameters from data.
%  [MU,KAPPA,RBAR,SD] = VM_ESTS(X,OPTIONS) estimates the angular mean (MU)
%  and concentration parameter (KAPPA) of the Von Mises distribution, the
%  mean vector length (RBAR), and the approximate standard deviation (SD)
%  of the normal distribution from data X.  X is assumed to be an n-by-2
%  matrix, where each row represents an x,y coordinate pair, which defines
%  a bearing from origin (0,0).  OPTIONS is an argument for FMINBND, for
%  instance as set by OPTIMSET.

n1 = size(x,1);
lx = sqrt(sum(x.^2,2));                % Get length of each vector
idx = find(lx);                        % Get rid of 0-length vectors
n = length(idx);
if n<n1,
  x = x(idx,:);
end  
x = x./repmat(lx,1,2);                 % Make unit vectors
r = sum(x);
mux = r/n;                             % Mean x,y coordinate
mu = 180/pi*atan2(mux(1),mux(2));      % Mean angle (use y/x for math convention)
Rbar = norm(mux);                      % Length of mean vector
sd = 180/pi*sqrt(2*(1-Rbar));          % Angular standard deviation
if Rbar<=(1/sqrt(n)),
  kappa = 0;
else
  kappa = fminbnd('diffkr',eps,5e5,options,Rbar,n);
end  
