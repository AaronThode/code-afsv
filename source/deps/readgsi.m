% readsgi
% matlab script to read and display sample of gsi file.
% whose name is in fn
%function [x,t,head]=readgsi(fn,ctstart,tlen,formatt);
%  fn-file name, including extension
%  ctstart: c-time of start time...
%  nsec: number of seconds to read in...
%  formatt: 'datenum' or 'ctime'

function [omi,t,head]=readgsi(fn,ctstart,tlen,formatt)

t=[];omi=[];
head=readgsif_header(fn);
if nargin==1
    return
end
if nargin<4
    formatt='ctime';
end

if strcmp(formatt,'datenum')&ctstart>0
    %tfile_start=datenum(1970,1,1,0,0,head.ctbc);
    ctstart=86400*(ctstart-datenum(1970,1,1,0,0,0));
elseif strcmp(formatt,'seconds')&ctstart>0
    ctstart=-ctstart;  %Set to zero
end
%fid = fopen(char(fn),'r','ieee-le');
fs=head.Fs*(1+head.tdrift/86400);
fid = fopen(char(fn),'r','ieee-be');

if head.ctbc>ctstart&ctstart>0
    %disp('Start time less than file start');
    fclose(fid);
    return;
end

if ctstart<=0
    %disp('cstart less than or equal to zero; interpret as seconds offset from file start');
    tsec=abs(ctstart);
else
    tsec=ctstart-head.ctbc;
end

fseek(fid,512+head.nc*floor(tsec*fs)*2,-1);

omi = fread(fid,[head.nc floor(head.Fs*tlen)],'uint16');
try
    fclose(fid);
catch
    disp('readgsi.m: couldn''t close file');
    fclose('all');
end
t=(1:size(omi,2))/head.Fs;

%N = 128;
%spectrogram(omi,hanning(N),N/2,N,head.Fs,'yaxis') % gram
%end
 
if 1==0
   fn='S508A0T20080819T072522.gsi';
   tstart=0;
   %tstart=datenum(2008,8,19,12,0,0);
   tsec=3;
   [x,t,head]=readgsi(fn,tstart,tsec,'datenum');
   for I=1:size(x,1),
       subplot(size(x,1),1,I);
       spectrogram(x(I,:),256,128,256,1000,'yaxis')
       
       
   end
   
end