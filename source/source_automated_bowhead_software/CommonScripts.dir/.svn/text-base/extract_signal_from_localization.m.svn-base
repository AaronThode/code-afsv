%function [y,start_ctime,t,head]=extract_signal_from_localization(localization,Istation, Icall,filename,start_str,nchan,buffer_time)
%% Give a station and index number, extract original time input into
%% morphological processor.  Useful for debugging or cross-matching
%% purposes...
%%  start_str: if 'debug' using the 'debug_time' field

function [y,start_ctime,t,head]=extract_signal_from_localization(localization,Istation, Icall,filename,start_str,nchan,buffer_time)

if ~exist('nchan')
    nchan=1;
end
if strcmp(start_str,'debug')
    start_ctime=(localization{Icall}.ctime_debug(Istation))-buffer_time;
    tlen=mean(localization{Icall}.duration_debug(Istation))+2*buffer_time;

else
    start_ctime=localization{Icall}.ctime_min(Istation)-buffer_time;
    tlen=mean(localization{Icall}.Totalduration(Istation))+2*buffer_time;

end
[y,t,head]=readfile(filename,start_ctime,tlen,nchan,'ctime','calibrate');

