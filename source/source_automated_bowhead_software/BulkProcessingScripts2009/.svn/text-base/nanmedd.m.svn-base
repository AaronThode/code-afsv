function y = nanmedd(x,dim)
%NANMEDD NaN protected median value across any dimension
%   NANMEDD(X,DIM) returns the median treating NaNs as missing values. 
%   For vectors, NANMEDD(X) is the median value of the non-NaN
%   elements in X.  For matrices, NANMEDIAN(X) is a row vector
%   containing the median value of each column, ignoring NaNs.
%
%   NANMEDD(X,DIM) takes the median along the dimension DIM of X.
%
%   See also NANMEAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   Combines the features of NANMEDIAN which takes the NAN protected
%   median across rows only, with MEDIAN which takes the median
%   across any dimension but is not NAN protected.
%   Chris Nations

if nargin==1, 
  dim = min(find(size(x)~=1)); 
  if isempty(dim), dim = 1; end
end
if isempty(x) % Check for empty input.
  y = NaN;
  return
end

siz = [size(x) ones(1,dim-ndims(x))];
m = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),m,prod(siz)/m);

x = sort(x); % NaNs are forced to the bottom of each column

% Replace NaNs with zeros.
nans = isnan(x);
x(nans) = 0;

n = m-sum(nans);
y = zeros(size(n));

% Odd columns
odd = find(rem(n,2)==1 & n>0);
idx =(n(odd)+1)/2 + (odd-1)*m;
y(odd) = x(idx);

% Even columns
even = find(rem(n,2)==0 & n>0);
idx1 = n(even)/2 + (even-1)*m;
idx2 = n(even)/2+1 + (even-1)*m;
y(even) = (x(idx1)+x(idx2))/2;

% All NaN columns
i = find(n==0);
y(i) = i + nan;

% Permute and reshape back
siz(dim) = 1;
y = ipermute(reshape(y,siz(perm)),perm);
