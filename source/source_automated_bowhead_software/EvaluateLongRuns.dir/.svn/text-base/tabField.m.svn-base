%function fs = TabField(s,n)
%
%  Where s is a string containing Tab-delimited fields,
%  this function will return the nth field, with trailing
%  and leading spaces removed.  If n is an array, then this
%  function returns a cell array, in which each cell is
%  the field indicated by the corresponding element in n.

function fs = TabField(s,n)
 %mbcharvector(s);
 %mbint(n);
  if isempty(n)
   fs = '';
  elseif length(n) == 1
   xx = find((s == char(9)));

    if n > (length(xx) + 1)
     fs = [];
    else
      if isempty(xx)
       fs = s;
      else
       xx = [0;xx(:);length(s) + 1];
       fs = s((xx(n) + 1):(xx(n + 1) - 1));
      end
    end
   fs = strtrim(fs);
  else
   fs = cell(size(n));
    for ii = 1:prod(size(n))
     fs{ii} = TabField(s,n(ii));
    end
  end
return