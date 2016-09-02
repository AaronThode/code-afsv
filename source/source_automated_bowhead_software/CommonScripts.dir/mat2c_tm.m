% tc = MAT2C_TM(tm)
%
%  Where tm is a 6-element row vector containing a date/time in MATLAB's
%  format, this function will return in tc a scalar containing the same
%  date/time in the standard C format, as the number of seconds past
%  Midnight, January 1, 1970.
%
%
% [tc,fs] = MAT2C_TM(tm)
%
%  If the second output argument, fs, if given, then only an integer value
%  will be returned in tc, with the fractional part of the seconds in fs.
%


function [tc,fs] = mat2c_tm(tm)
tm=datevec(tm);
a = version;
 if (a(1) == 4) & (a(2) == '.')
  error('Missing Mat2C_TM.mex');
 else
  tc = datenum(tm(1),tm(2),tm(3),tm(4),tm(5),tm(6))*(24*60*60) - 6.216730560000000e+010;
   if nargout > 1
    fs = tc - fix(tc);
    tc = fix(tc);
   end
 end
