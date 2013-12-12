% rdsiof
% matlab script to read and display sample of sio file.
% whose name is in fn
%function [x,t,head,formatt]=readsiof(fn,ctstart,tlen);
%  fn-file name, including extension
%  ctstart: c-time of start time...
%  nsec: number of seconds to read in...

function [omi,t,head]=readsiof(fn,ctstart,tlen,formatt),

t=[];omi=[];
head=readsiof_header(fn);
if nargin==1,
    return
end
if nargin<4,
    formatt='ctime';
end

if strcmp(formatt,'datenum'),
    %tfile_start=datenum(1970,1,1,0,0,head.ctbc);
    ctstart=86400*(ctstart-datenum(1970,1,1,0,0,0));

end
%fid = fopen(char(fn),'r','ieee-le');
fs=head.Fs*(1+head.tdrift/86400);
fid = fopen(char(fn),'r','ieee-be');

if head.ctbc>ctstart,
    disp('Start time less than file start');
    return;
end

%if strcmp(formatt,'datenum');
%    tvec=datevec(tstart-tfile_start);
%    tsec=tvec(:,6)+tvec(:,5)*60+tvec(:,4)*3600+tvec(:,3)*24*3600;
%else
    tsec=ctstart-head.ctbc;
%end
fseek(fid,512+floor(tsec*fs)*2,-1);

% read first 60,00 samples of time series (1 min)
omi = fread(fid,floor(head.Fs*tlen),'uint16');
fclose(fid);

t=(1:length(omi))/head.Fs;

%N = 128;
%spectrogram(omi,hanning(N),N/2,N,head.Fs,'yaxis') % gram
%end
 
if 1==0,
   fn='D06_1B1_20060921T000000.sio';
   tstart=0;
   tsec=60;
   [x,head]=readsiof(fn,tstart,tsec);
end