%%%%read_Java_SEL.m%%%%%
%  Read a *.sel file produced by PreProcess_PSD_detection java program.
% function [SEL,Tsec,Tabs,params]=read_Java_SEL(fname,tabs_start,tduration)
% Input:
%
% fname-filename including absolute path and .sel extension.
%       Example formats:
%       Single unbeamed channel: (only the part from "tstart" String matters):
%%%     pm_Sound_2007_07210719_tstart2007_07_21_07_19_00_Fs4096_100to1000Hz_Nfft1024_dn512.sel
%
%       Beamed channel (newer version; contains header..):
%           10_tstart2009_08_22_11_08_56_100to300Hz_5beams.sel
%
% tabs_start.  If greater than zero, interpret as a datenumber.
%   If less than or equal to zero, then read data in starting at abs(tabs_start)
%   seconds into file.
%
% tduration:  Time duration of data desired in seconds
%
%  Output:
%
%  SEL(Nbeam,Ntime): sound exposure level in linear units.  May be calibrated or
%  uncalibrated depending on original source file
%
%  Tsec:  Associated time points in seconds from start of file.  The time
%  origin is actually Nfft/Fs sec from the file start contained in the
%  fname
%
%  Tabs: If absolute start time is recorded, contains datenumbers of
%  appropriate SEL..
%
%  params: information used to compute SEL, including Nfft, etc:
%       params.fmin=fmin;
%       params.fmax=fmax;
%       params.Nfft=Nfft;
%       params.Fs=Fs;
%       params.dn=dn;
%       params.tstart_file;
%       params.nchan
%  params.angles=angles;


%Updated Jan 14, 2010 to include beamformed SEL file style, which will have
%       'beams' in the file name

function [SEL,Tsec,Tabs,params]=read_Java_SEL(fname,tabs_start,tduration)

%%Parse filename
params.angles=[];
bytes_sample=4;  %bytes per sample. A float is 4 bytes...
if isempty(findstr(fname,'beams'))
    [Fs,fmin,fmax,Nfft,dn,ichan]=load_single_channel_params(fname);
    nchan=1;
    byte_offset_header=0;
else
   [Fs,fmin,fmax,Nfft,dn,angles]=load_beamed_channel_params(fname);
   nchan=length(angles);
   byte_offset_header=512;
   params.angles=angles;
   ichan=0;
end
tstart_file=getTstartFileName(fname);
fid=fopen(fname,'r','ieee-be');
status=fseek(fid,byte_offset_header,-1);



dT=dn/Fs;  %Time increment for energy function
%Time of first Sel point actually assigned at midpoint of first Nfft

tstart_file=tstart_file+datenum(0,0,0,0,0,0.5*Nfft/Fs);
samples_read=floor(tduration/dT);  %Number of energy samples read...
%%Determine byte offset if want to load at non-zero start
if tabs_start<=0, % If we desire offset in seconds from file start..
    %disp('tabs_start interpreted as seconds offset');
    byte_offset=nchan*bytes_sample*floor(abs(tabs_start)/dT);
elseif tstart_file>tabs_start
    error(sprintf('Time requested: %s is less than file start time %s', ...
        datestr(tabs_start),datestr(tstart_file)));
else
    %disp(' date information available');
    tabs_offset=datevec(tabs_start-tstart_file);
    tabs_offset=(24*3600*tabs_offset(:,3)+3600*tabs_offset(:,4)+60*tabs_offset(:,5)+tabs_offset(:,6));
    byte_offset=nchan*bytes_sample*floor(tabs_offset/dT);
end

status=fseek(fid,byte_offset,0);
%keyboard;
if status==0
    SEL=fread(fid,[nchan samples_read],'float32'); 
else
    disp('Requested time past end of file');
    SEL=[];
end
if isempty(SEL)
    disp('EOF reached, result dumped..');
end
%Create time axis..assign center off FFT sample as representative time...
Tsec=dT/2+(byte_offset/(bytes_sample*nchan))*dT+(0:(size(SEL,2)-1))*dT;
if tabs_start<0,
    %The Nfft/2 reflects fact that energy SEL
    %sample uses Nfft data points, so place time of estimate halfway
    %into data sample
    Tabs=[];
else
    Tabs=tstart_file+datenum(0,0,0,0,0,Tsec);
    
end
fclose(fid);

%%Clear out trailing zeros in case of 'Inf' time sample request
Igood=find(SEL>0);
SEL=SEL(Igood);
Tabs=Tabs(Igood);
Tsec=Tsec(Igood);

params.fmin=fmin;
params.fmax=fmax;
params.Nfft=Nfft;
params.Fs=Fs;
params.dn=dn;
params.tstart_file=tstart_file;
params.nchan=nchan;
params.ichan=ichan;
%params.angles=angles;



end

function tstart_file=getTstartFileName(fname)
Is=findstr(fname,'tstart')+6;
year=str2double(fname(Is+(0:3)));
mon=str2double(fname(Is+(5:6)));
day=str2double(fname(Is+(8:9)));
hr=str2double(fname(Is+(11:12)));
min=str2double(fname(Is+(14:15)));
sec=str2double(fname(Is+(17:18)));
tstart_file=datenum(year,mon,day,hr,min,sec);
end

function [Fs,fmin,fmax,Nfft,dn,angles]=load_beamed_channel_params(fname)
fid = fopen(char(fname),'r','ieee-be'); %for s1sT
Fs = fread(fid,1,'double');
fmin = fread(fid,1,'int32');
fmax = fread(fid,1,'int32');
Nfft = fread(fid,1,'int32');
dn = fread(fid,1,'int32');
%nstart=fread(fid,1,'int32');
nchan=fread(fid,1,'int32');
angles=fread(fid,nchan,'double');
fclose(fid);

end

function [Fs,fmin,fmax,Nfft,dn,ichan]=load_single_channel_params(fname)


Is=findstr(fname,'Fs')+2;
Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
Fs=str2double(fname(Is:Ie));

Is=Ie+2;
Ie=findstr(fname(Is:end),'to')-2;Ie=Is+Ie(1);
fmin=str2double(fname(Is:Ie));

Is=Ie+3;
Ie=findstr(fname(Is:end),'Hz')-2;Ie=Is+Ie(1);
fmax=str2double(fname(Is:Ie));

Is=findstr(fname,'Nfft')+4;
Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
Nfft=str2double(fname(Is:Ie));

Is=findstr(fname,'dn')+2;
Ie=findstr(fname(Is:end),'.sel')-2;Ie=Is+Ie(1);
dn=str2double(fname(Is:Ie));

%Delphine commented out--may have to uncomment 
% Is=findstr(fname,'chann')+5;
% Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
% ichan=str2double(fname(Is:Ie));
ichan=1;
end


