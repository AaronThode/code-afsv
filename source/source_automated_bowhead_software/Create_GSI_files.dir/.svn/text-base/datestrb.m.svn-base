function date = datestrs(when,prec,fullyear)
%  Where when is a 6-element array where the elements of when represent the
%  year, month, day, hour, minute, and second, in that order, datestrb(when)
%  will return a string which represents the date and time.
%
%  datestrb(when,prec), where prec is a positive integer, will return a string
%  which shows the seconds to prec places after the decimal point.

 if nargin < 3
  fullyear = 0;
 end

 if length(fullyear) ~= 1
  fullyear = 0;
 end

 if isnan(fullyear);
  fullyear = 0;
 end

MONTHLIST = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
fracSec = when(6) - fix(when(6));
when = fix(when);

s = num2str(when(3));
 if size(s,2) == 1
  s = ['0' s];
 end
date = [s ' ' MONTHLIST(when(2),:) ' '];

 if fullyear
  s = num2str(rem(when(1),10000));
   while size(s,2) < 4
    s = ['0' s];
   end
  date = [date s '  '];
 else
  s = num2str(rem(when(1),100));
   while size(s,2) < 2
    s = ['0' s];
   end
  date = [date s '  '];
 end

s = num2str(when(4));
 if size(s,2) == 1
  s = ['0' s];
 end
date = [date s ':'];

s = num2str(when(5));
 if size(s,2) == 1
  s = ['0' s];
 end
date = [date s ':'];

s = num2str(when(6));
 if size(s,2) == 1
  s = ['0' s];
 end
date = [date s];

 if nargin >= 2
   if isstr(prec)
    prec = fix(str2num(prec));
   end

   if prec >= 1
    s = sprintf('%.*f',prec,fracSec);
     while(s(1) ~= '.');
      s = s(2:length(s));
     end
    date = [date s];
   end
 end
