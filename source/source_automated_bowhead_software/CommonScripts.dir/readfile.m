%function [y,t,head]=readfile(rawfile,cbegin,tlen,nchan,formatt,calibrate)
% Input Parameters:
% rawfile = Include extension;
% cbegin = Start time;
% tlen = Length of sample to load;
% nchan = Index of desired channel, 1 for sound;
% formatt = Describes time input, string 'ctime' or 'datenum'; 
% calibrate = String 'calibrate' to convert Volts to microP;

function [y,t,head]=readfile(rawfile,cbegin,tlen,nchan,formatt,calibrate)
if findstr(rawfile,'.sio')

    [y,t,head]=readsiof(rawfile,cbegin,tlen,formatt);
    y=y-0.5*(2^16);  %Remove DC bias in A/D converter

elseif findstr(rawfile,'.gsi')
    [y,t,head]=readgsi(rawfile,cbegin,tlen,formatt);
    if isempty(y)  %request time is befine file start
        dt=cbegin-head.ctbc;
        tlen=tlen+dt;
        [y,t,head]=readgsi(rawfile,0,tlen,formatt);
    
    end
    
        
    y=y(nchan,:);
    y=calibrate_GSI_signal(y, 'DASARC');
        
    
    if (abs(max(size(y))-floor(head.Fs*tlen))>2),
        disp('End of file reached, setting y to empty');
        y=[];
    end
end
