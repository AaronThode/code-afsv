% tm = c2mat_tm(tc)
%
%  Where tc is a scalar containing a date/time as stored by C, as the
%  number of seconds past Midnight, January 1, 1970, this function will
%  return in tm a six-element row vector containing the same date/time
%  in MATLAB's format.
%
%
% tm = c2mat_tm(tc,fs)
%
%  The additional argument, fs, will be added to the seconds field of the
%  result.  This is to allow times to be specifed to fractions of a second.
%  If no second input argument is given, then the fractional part of tc, if
%  present, will be used.
%

function tm = c2mat_tm(tc,fs)
a = version;
 if (a(1) == 4) & (a(2) == '.')
  error('Missing C2Mat_TM.mex');
 else
   if nargin > 2
    tc = tc + fs;
   end
  tm = datevec((tc+6.216730560000000e+010)/(24*60*60));
 end