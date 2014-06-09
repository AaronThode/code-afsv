%%%%read_Java_PSD.m%%%%%
%  Read a *.psd file produced by PreProcess_PSD_detection java program.
% function [PSD,F,Tsec,Tabs,params]=read_Java_PSD(fname,tabs_start,tduration,header_only)
% Input:
%
% fname-filename including absolute path and .psd extension.
%       Example formats:
%       Single unbeamed channel: (only the part from "tstart" String matters):
%   %%%pm_Sound_2007_07210719_tstart2007_07_21_07_19_00_chann0_Fs4096_Nfft1024_ovlap512_Nsamps39.psd
%
%       Beamed channel (newer version; contains header..):
%           10_tstart2009_08_22_11_08_56_5beams.psd
%
% tabs_start.  If greater than zero, interpret as a datenumber.
%   If less than or equal to zero, then read data in starting at abs(tabs_start)
%   seconds into file.
%
% tduration:  Time duration of data desired in seconds.  Can be Inf
%
%   header_only:  if exists, return params only and emptys for other
%   variables
%  Output:
%
%  PSD(Nfreq,Nbeam,Ntime): two-sided pressure spectral density in linear units.  May be calibrated or
%  uncalibrated depending on original source file.  Lowest frequency (DC)
%  is first row, up to Nfft/2 rows
%
%  F: associated frequency bins...
%
%  Tsec:  Associated time points in seconds from start of file.  The time
%  origin is actually Nfft/Fs sec from the file start contained in the
%  fname
%
%  Tabs: If absolute start time is recorded, contains datenumbers of
%  appropriate PSD..
%
%  params: information used to compute PSD, including Nfft, etc:
%  params.Nfft=Nfft;
%  params.Fs=Fs;
%  params.dn=dn;    %samples advanced between FFT computations
%  params.Nsamps=Nsamps;  %Averaged samples
%  params.Nmax...     %%Number of frequency bins loaded (may be less than Nfft/2)
%  params.tstart_file=tstart_file;
%  params.nchan=nchan;
%  params.angles=angles;
%  params.total_time_sec:  %Total time covered by PSD file

%Updated Jan 14, 2010 to include beamformed PSD file style, which will have
%       'beams' in the file name

function [PSD,F,Tsec,Tabs,params]=read_Java_PSD(fname,tabs_start,tduration,header_only)

if nargin<4
    header_only=0;
end
%%Parse filename
params.angles=[];
bytes_sample=4;  %bytes per sample. A float is 4 bytes...
if isempty(findstr(fname,'beam'))
    [Fs,Nfft,dn,Nsamps,ichan,Nmax]=load_single_channel_params(fname);
    nchan=1;
    byte_offset_header=0;
    beam_flag=0;
else
    [Fs,Nfft,dn,Nsamps,angles]=load_beamed_channel_params(fname);
    nchan=length(angles);
    byte_offset_header=8*512;
    params.angles=angles;
    %Nangles=length(angles);
    beam_flag=1;
    ichan=0;
end
tstart_file=getTstartFileName(fname);
fid=fopen(fname,'r','ieee-be');

%%Added Oct. 2012:  Check size of file and return number of samples.
fseek(fid,0,1);  %Advance to end of file
nbytes=ftell(fid);
fseek(fid,0,-1);  %rewind the file

    
%%If beamformed data, advance past header file
status=fseek(fid,byte_offset_header,-1);

%%Error in function below, dT should be actually  dn*(Nsamps-1) only
%%dT=(Nfft+(dn-1)*Nsamps)/Fs;  %Time increment for PSD function

dT=dn*Nsamps/Fs;
Tstart=dT/2;  %State that time assigned to a PSD is mid-point of averaging time.
%Npsd=1+Nfft/2;  %Number of points to read from Psd, plus one check point
Npsd=1+Nmax;

%Time of first PSD point actually assigned at midpoint of sample averages
tstart_file=tstart_file+datenum(0,0,0,0,0,0.5*dT);
samples_read=floor(tduration/dT);  %Number of energy samples read...

params.total_time_sec=dT*nbytes/(nchan*Npsd*bytes_sample);
params.tend_file=tstart_file+datenum(0,0,0,0,0,params.total_time_sec);
params.Nfft=Nfft;
params.Fs=Fs;
params.dn=dn;
params.Nsamps=Nsamps;
params.tstart_file=tstart_file;
params.nchan=nchan;
params.ichan=ichan;
params.Nmax=Nmax;

if header_only~=0
    PSD=[];F=[];Tsec=[];Tabs=[];
    fclose(fid);

    return
end

%%Determine byte offset if want to load at non-zero start
if tabs_start<=0, % If we desire offset in seconds from file start..
    %disp('tabs_start interpreted as seconds offset');
    byte_offset=nchan*Npsd*bytes_sample*floor(abs(tabs_start)/dT);
    
elseif tstart_file>tabs_start
    error(sprintf('Time requested: %s is less than file start time %s', ...
        datestr(tabs_start),datestr(tstart_file)));
else
    %disp(' date information available');
    tabs_offset=datevec(tabs_start-tstart_file);
    tabs_offset=(24*3600*tabs_offset(:,3)+3600*tabs_offset(:,4)+60*tabs_offset(:,5)+tabs_offset(:,6));
    byte_offset=nchan*bytes_sample*Npsd*floor(tabs_offset/dT);
end

status=fseek(fid,byte_offset,0);
tmp=zeros(Npsd,1);
%keyboard;
if status==0
    PSD=fread(fid,[Npsd samples_read],'float32');
    if ndims(PSD)==3
        PSD=PSD(2:end,:,:);
    else
        PSD=PSD(2:end,:);
    end
else
    disp('Requested time past end of file');
    PSD=[];
end
if isempty(PSD),
    disp('EOF reached, result dumped..');
end
%Create time axis..assign center off FFT sample as representative time...
F=(0:(Nfft/2))*(Fs/Nfft);
try
    F=F(1:Nmax);
catch
    disp('F may be inaccurate, autocorrelation output perhaps?')
end
Tsec=Tstart+(byte_offset/(Npsd*bytes_sample*nchan))*dT+(0:(size(PSD,2)-1))*dT;

% if tabs_start<0
%     %The Nfft/2 reflects fact that energy SEL
%     %sample uses Nfft data points, so place time of estimate halfway
%     %into data sample
%     Tabs=tstart_file
% else
    Tabs=tstart_file+datenum(0,0,0,0,0,Tsec);
    
% end
fclose(fid);

try
    params.angles=angles;
end

end

%Retreive function name from input string
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

function [Fs,Nfft,dn,Nsamps,angles]=load_beamed_channel_params(fname)
fid = fopen(char(fname),'r','ieee-be'); %for s1sT
Fs = fread(fid,1,'double');
Nfft = fread(fid,1,'int32');
dn = fread(fid,1,'int32');
Nsamps = fread(fid,1,'int32');
%nstart=fread(fid,1,'int32');
%nchan=fread(fid,1,'int32');
%angles=fread(fid,nchan,'double');
angles=fread(fid,1,'double');
fclose(fid);

end

function [Fs,Nfft,dn,Nsamps,ichan,Nmax]=load_single_channel_params(fname)


Is=findstr(fname,'Fs')+2;
Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
Fs=str2double(fname(Is:Ie));

Is=findstr(fname,'Nfft')+4;
Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
Nfft=str2double(fname(Is:Ie));

Is=findstr(fname,'dn')+2;
Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
dn=str2double(fname(Is:Ie));
%dn=Nfft-ovlap;

Is=findstr(fname,'chann')+5;
Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
ichan=str2double(fname(Is:Ie));


Is=findstr(fname,'Nmax')+4;
if ~isempty(Is)
    Ie=findstr(fname(Is:end),'.psd')-1;Ie=Is+Ie(1);
    Nmax=str2num(fname(Is:Ie));
    
    Is=findstr(fname,'Nsamps')+6;
    Ie=findstr(fname(Is:end),'_')-2;Ie=Is+Ie(1);
    Nsamps=str2double(fname(Is:Ie));
    
else
    Nmax=Nfft/2;
    
    Is=findstr(fname,'Nsamps')+6;
    Ie=findstr(fname(Is:end),'.psd')-2;Ie=Is+Ie(1);
    Nsamps=str2double(fname(Is:Ie));
    
end
end
