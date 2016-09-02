%function [y,start_ctime,t,head]=extract_signal_from_station(station,Idetect,filename,start_str,nchan)
%% Give a station and index number, extract original time input into
%% morphological processor.  Useful for debugging or cross-matching
%% purposes...
%%  start_str: if 'debug' using the 'debug_time' field

function [y,start_ctime,tlen,head]=extract_signal_from_station(station,Idetect,filename,start_str,nchan)

if ~exist('nchan')
    nchan=1;
end
if strcmp(start_str,'debug')
    start_ctime=(station.ctime_debug(Idetect));
    tlen=mean(station.duration_debug(Idetect));

else
    start_ctime=station.ctime_min(Idetect);
    tlen=mean(station.Totalduration(Idetect));

end
[y,t,head]=readfile(filename,start_ctime,tlen,nchan,'ctime','calibrate');

