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
for i = 1:B
    %U = ceil(n .* rand(n,1));          % The guts of UNIDRND.M without error checking.
    %xb = x(U,:);
    mux_hat(i,:) = vm_ests_uneq(x(U(:,i),:),options,flag); % Estimation accounting for lengths.
end

[~,kappa,sd] = vm_ests(mux_hat,options);     % Estimation ignoring lengths.
%toc
end
