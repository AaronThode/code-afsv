function [y,t,head]=readGSIfile(rawfile,cbegin,tlen,nchan,formatt,calibrate)
%function [y,t,head]=readGSIfile(rawfile,cbegin,tlen,nchan,formatt,calibrate)
% Input Parameters:
% rawfile = Include extension;
% cbegin = Start time;
% tlen = Length of sample to load;
% nchan = Index of desired channel, 1 for sound;
% formatt = Describes time input, string 'ctime' or 'datenum';
% calibrate = String 'calibrate' to convert Volts to microP;

if strfind(rawfile,'.sio')
    
    [y,t,head]=readsiof(rawfile,cbegin,tlen,formatt);
    y=y-0.5*(2^16);  %Remove DC bias in A/D converter
    
elseif strfind(rawfile,'.gsi')
    [y,t,head]=readgsi(rawfile,cbegin,tlen,formatt);
    if isempty(y)  %request time is befine file start
        dt=cbegin-head.ctbc;
        tlen=tlen+dt;
        [y,t,head]=readgsi(rawfile,0,tlen,formatt);
        
    end
    y=y(nchan,:);
    if strcmp(calibrate,'calibrate')
        y=y-0.5*(2^16);  %Remove DC bias in A/D converter
        
        y=y*(2.5/65535)*(10^(150/20));
    end
    y=y.';
    
    if (abs(size(y,1)-floor(head.Fs*tlen))>2)
        disp('Warning: end of file reached before requested data acquired');
        %y=[];
    end
end

end  %function readGSIfile
