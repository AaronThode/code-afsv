
function [x,t,head,tdate]=load_mt(fname,tstart,tsec,uncal_chc)
%function [x,t,head,tdate]=load_mt(fname,tstart,tsec,uncal_chc),
% fname: file name
% tstart: start of file in seconds
% tsec: length of segment desired in seconds
% uncal_chc; If exists, return uncalibrated values.
% All times in seconds

if isempty(strfind(lower(fname),'.mt')),
    fname=[fname '.mt'];    
end
head=read_mt_header(fname);
Fs=head.Fs;
if strcmp(head.signing,'S'),
    dattype=['int' num2str(8*head.wordsize)];
else
    dattype=['uint' num2str(8*head.wordsize)];    
end
if strcmp(head.swapping,'U'),
    fid=fopen(fname,'r','ieee-be');
else
    fid=fopen(fname,'r','ieee-le');
end

if fid<0,
    
    keyboard;
end
fseek(fid,512,-1);  %Skip the 512 bytes

offset=head.wordsize*floor(Fs*tstart);
status=fseek(fid,offset,0);
x=fread(fid,floor(tsec*Fs)+1,dattype);
fclose(fid);
%Optional calibration...
N=(2.^head.samplebits)-1;
if ~exist('uncal_chc'),
    switch head.signing
        case 'S',
            calmean=0.5*(head.calmin+head.calmax);
            x=calmean+x*(head.calmax-head.calmin)/N;
        case 'U',
            calmean=head.calmin;
            x=calmean+x*(head.calmax-head.calmin)/N;
    end
    if findstr(head.abbrev,'Sound')
        x=x*1000; %convert from mPa to microPa
    end
end
t=(0:(length(x)-1))/Fs;
tdate=datenum(head.year,head.month,head.day,head.hours,head.minutes,tstart+t+head.seconds+head.msec/1000);
