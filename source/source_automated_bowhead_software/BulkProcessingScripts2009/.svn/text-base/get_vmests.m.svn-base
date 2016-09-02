function [mu,kappa,sd] = get_vmests(x,B)
%GET_VMESTS Calculates bootstrap estimate of mean and standard error of bearings.
%  [MU,KAPPA,SD] = GET_VMESTS(X,B) assumes that X is an n-by-2 matrix of x-y
%  coordinate pairs and B is a scalar denoting the number of bootstrap iterations.
%  MU is the mean bearing, KAPPA is an estimate of the standard error of the mean
%  expressed as the von Mises concentration parameter, and SD is the standard error
%  estimate in degrees expressed as for a linear (not a circular variable).
%
%  Note: B may be reduced for greater speed.

n = size(x,1);
flag = 0;           % Set to 1 to calculate estimates of kappa and Rbar for each
%   bootstrap sample, 0 for mean vector only.
mux_hat = zeros(B,2);
options = optimset('Display','off');
mux = vm_ests_uneq(x,options,flag);             % Get mean angle from the data.
mu = 180/pi*atan2(mux(1),mux(2));
%tic
%for i = 1:B,
%  U = ceil(n .* rand(n,1));          % The guts of UNIDRND.M without error checking.
%  xb = x(U,:);
%  mux_hat(i,:) = vm_ests_uneq(xb,options,flag); % Estimation accounting for lengths.
%end

U = ceil(n .* rand(n,B));          % The guts of UNIDRND.M without error checking.
for i = 1:B,
    %U = ceil(n .* rand(n,1));          % The guts of UNIDRND.M without error checking.
    %xb = x(U,:);
    mux_hat(i,:) = vm_ests_uneq(x(U(:,i),:),options,flag); % Estimation accounting for lengths.
end

[junk,kappa,sd] = vm_ests(mux_hat,options);     % Estimation ignoring lengths.
%toc
end

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

%n1 = size(x,1);

lx = sqrt(sum(x.^2,2));                % Get length of each vector
idx = find(lx);                        % Get rid of 0-length vectors
%n = length(idx);
%if n<n1,
x = x(idx,:);
%end
r = sum(x);
mux = r/sum(lx);                       % Mean x,y coordinate


end


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
end

% %tic
% U=ceil(n.*rand(n,B));
% for I=1:B
%    mux_hat(I,:)=vm_ests_uneq(x(U(:,I)),options,flag);
% end
% [junk,kappa,sd] = vm_ests(mux_hat,options);     % Estimation ignoring lengths.
% %toc

