%%%%%%function ctime_out=datenum2ctime(twant,ctime_in)%%%
%% ctime_in: Ctime of a date within 24 hours of desired ctimes
%% twant: datenumber of desired output ctime

function ctime_out=datenum2ctime(twant,ctime_in)

tabs_in=datenum(1970,1,1,0,1,ctime_in);
tdiff=datevec(tabs_in-twant);
tdiff=tdiff(6)+60*tdiff(5)+3600*tdiff(4)+24*3600*tdiff(3);
ctime_out=ctime_in-tdiff;