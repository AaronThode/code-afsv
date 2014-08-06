function c = cond2(A)
%COND2   Condition number with respect to inversion.
%   COND2(X) returns the 2-norm condition number (the ratio of the
%   largest singular value of X to the smallest).  Large condition
%   numbers indicate a nearly singular matrix.
%
%   Modified to handle only the 2-norm, eliminate some error checking,
%   and other condition testing.
%      Chris Nations  8/20/02

s = svd(A);
if any(s == 0)   % Handle singular matrix
  c = Inf;
else
  c = max(s)./min(s);
end
