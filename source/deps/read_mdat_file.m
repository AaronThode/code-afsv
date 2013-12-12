%%%%%%read_mdat_file.m%%%%%%%%%%%
% [x,fparms]=read_mdat_file(fname_in,tstart,nsec,calibrate_acoustic)
% Read in a MDAT file.  WARNING!  LOG file must be in same location.
%
% Input:
% fname: name of *MDAT file, including full pathname.  MUST include 'MDAT' extension!
% tstart-datenumber of desired start time.  If zero, data are read from
%   start of file
% nsec-number of seconds to be pulled.
% calibrate_acoustic: if 1, then output in terms of acoustics, if 0 output in terms of voltage
%       if 2, output raw count.
%       Default is raw count
% 
% Output:
%   x:  [Ntimes Nsensor] output array
% Output of fparms:
%   tfs: datenumber of the start of the file
%   tfe: datenumber of end of file
%   fs: sampling rate in Hz
%   nc: number of channels
%   sens:  sensitivity in V/uPa
%   Igood: 1 if channel is good, 0 if bad
%   nstart:  number of time samples into MDAT collected 
%   nsamples:   number of time samples read
%   synch: structure of time-synchronization with other files...
%       file (string of relative pathname to synched file)
%       time (datenum of synchronization time)
%       offset (seconds, measured at synch time)
%       drift (msec/hour) measured relative to synch time
%   geom:  structure of geometry:
%       D:  water depth in m
%       rd: channel depths in meters
%       spacing: channel separation in meters from next (higher channel)
%       tiltx, tilty: tilt of array in degrees

% Version 3, Aug 2013--loads in time synch and geometry from LOG file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Example of Required format for a MDAT LOG file: (synch optional):
%
%
% Acoustic Data Logger deployment settings
% Run time: 7200
% Sleep time: 7200
% Channels: 8
% Sample rate: 6250 Hz per channel
% 
% Start time = Tue Aug 31 10:00:20 2010
%   0.802225 seconds
% Tilt-X:  13.8   Tilt-Y:  -6.1  degrees
% Stop time = Tue Aug 31 12:00:20 2010
%   0.002700 seconds
% Battery Voltage:  9.84
% 
% Channel depths (m): 50.9106   47.9134   44.8654   41.8301   38.8583   35.9627   32.8385   30.1207
% Channel spacing (m): 0         2.9972    3.0480    3.0353    2.9718    2.8956    3.1242    2.7178
% Channel quality: 1     1     1     1     1     1    1     1     
% Water depth: 55 m
% 
% Synch file: ../Top_Unit/103.MDAT
% Synch time = Tue Aug 31 10:00:29 2010
%   0.0 seconds
% Time offset: 2.50596 seconds
% Time drift:  9.893566322478158 msec/hour


function [x,fparms]=read_mdat_file(fname_in,tstart,ns,calibrate_acoustic)

[dirname,fname,extt] = fileparts(fname_in);
fname=[fname extt];
dirname=dirname(~isspace(dirname));  %deblank only works on trailing edges
if isempty(dirname)
    dirname='.';
end

if ~exist('calibrate_acoustic')
    calibrate_acoustic=2;
end
%fs=50000;
VREF=2.5;
Nmax=(2^16)-1;
bias=(Nmax+1)/2;
if calibrate_acoustic<2
    scale=(0.038/1.00)*VREF/Nmax;  %An input voltage of 0.05 V p-p produces a 1.2 V p-p at A/D
elseif calibrate_acoustic==2
    scale=1;
end

%nc=1;
%np=ns*fs;

%Get data from log file
mydir=pwd;
cd(dirname)
[tfs,tfe,fs,nc,sens,Igood,geom,synch]=load_mdat_header(fname);

fparms=struct('tfs',tfs,'tfe',tfe,'fs',fs,'nc',nc,'geom',geom,'synch',synch,'sensitivity',sens,'Igood',Igood);
fid=fopen(fname,'r','ieee-be');
cd(mydir);

try
    if tstart>tfs&&tstart<(tfe-datenum(0,0,0,0,0,ns))
        offset=datevec(tstart-tfs);
        offset=(offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file-DO NOT ROUND!
        disp(sprintf('%i second offset from file start',offset));
    else
        %error('cannot read this time from this file');
        disp('desired start is before file start, setting offset to 0');
        offset=0;
        
        %offset=4;
    end
catch
    disp('File does not contain wanted time, will interpret tstart as elasped time in seconds');
    %keyboard;
    offset=tstart;
end

fparms.nstart=round(offset*fs);
%scale=1;
fskip=2*nc*fparms.nstart;
np=round(fs*ns);
fparms.nsamples=np;

res=fseek(fid,fskip,-1);
if res<0
    keyboard;
end
disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));
x0=fread(fid,[nc np],'uint16');
fclose(fid);
x=scale*(x0-bias);


%Crude acoustic level calibration
 if calibrate_acoustic==1
     disp('Acoustic calibration executed');
     %sens=182; %dB re 1uPa/V
     x=x.*(10^(-sens/20));
 end

%Future calibration corrections
x=x.';
cd(mydir);

end

function [tfs,tfe,fs,nc,sensitivity,Igood,geom,synch]=load_mdat_header(fname)
%[tfs,tfe,fs,nc,sens,Igood,geom,synch]=load_mdat_header(dirname,fname);
%   geom:  structure of geometry:
%       D:  water depth in m
%       rd: channel depths in meters
%       spacing: channel separation in meters from next (higher channel)
%       tiltx, tilty: tilt of array in degrees
% synch parameters:
%       file (string of relative pathname to coordinated data)
%       time (datenum of synchronization time)
%       offset (seconds, measured at synch time)
%       drift (msec/hour) measured relative to synch time
%       

fpre=findstr(lower(fname),'.mdat')-1;
log_name=[fname(1:fpre) '.LOG'];
fid=fopen(log_name,'r');
sensitivity=-160; %dB re 1uPa/V, default sensitivity
Igood=[];
synch=struct('file',[],'time',[],'offset',[],'drift',[]);
geom=struct('rd',[],'D',[],'spacing',[],'tiltx',[],'tilty',[]);

while 1
    tline=fgetl(fid);
    if ~ischar(tline),break,end
    tline_org=tline;
    tline=lower(tline);
    if findstr(tline,'channels:')>0,
        Idot=1+findstr(tline,':');
        nc=str2num(tline(Idot:end));
    elseif findstr(tline,'sample rate:')>0
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'hz')-1;
        fs=str2num(tline(Idot:Iend));
    elseif findstr(tline,'sensitivity:')>0
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'db')-1;
        sensitivity=str2num(tline(Idot:Iend));
    elseif findstr(tline,'start time')>0,
        Idot=findstr(tline,'=');
        tfs=parse_date(tline(Idot:end));
        tline=fgetl(fid);
        Iend=findstr(tline,'seconds')-1;
        tfs=tfs+datenum(0,0,0,0,0,str2num(tline(1:Iend)));
    elseif findstr(tline,'stop time')>0
        Idot=findstr(tline,'=');
        tfe=parse_date(tline(Idot:end));
        tline=fgetl(fid);
        Iend=findstr(tline,'seconds')-1;
        tfe=tfe+datenum(0,0,0,0,0,str2num(tline(1:Iend)));
    elseif findstr(tline,'synch time')>0
        Idot=findstr(tline,'=');
        synch.time=parse_date(tline(Idot:end));
        tline=fgetl(fid);
        Iend=findstr(tline,'seconds')-1;
        synch.time=synch.time+datenum(0,0,0,0,0,str2num(tline(1:Iend)));
    elseif findstr(tline,'tilt-x')>0,
        Idot=findstr(tline,':')+1;
        Iend=findstr(tline,'tilt-y')-1;
        geom.tiltx=str2num(tline(Idot:Iend));
        Idot=Iend+8;
        Iend=findstr(tline,'degrees')-1;
        geom.tilty=str2num(tline(Idot:Iend));
    elseif findstr(tline,'synch file:')>0
        Idot=1+findstr(tline,':');
        synch.file=(tline_org(Idot:end));
    elseif findstr(tline,'time offset:')>0
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'sec')-1;
        synch.offset=str2num(tline(Idot:Iend));
    elseif findstr(tline,'time drift:')>0
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'ms')-1;
        synch.drift=str2num(tline(Idot:Iend));
    elseif findstr(tline,'channel depths')>0
        Idot=1+findstr(tline,':');
        geom.rd=str2num(tline(Idot:end));
    elseif findstr(tline,'channel spacing')>0
        Idot=1+findstr(tline,':');
        geom.spacing=str2num(tline(Idot:end));
    elseif findstr(tline,'water depth:')>0
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'m')-1;
        geom.D=str2num(tline(Idot:Iend));
    elseif findstr(tline,'channel quality:')>0
        Idot=1+findstr(tline,':');
        Igood=str2num(tline(Idot:end));
    end
    
end

fclose(fid);
    function tabs=parse_date(str)
        
        Istart=findstr(str,'=');
        year=str2num(str((end-4):end));
        tm=datenum(str((end-13):(end-5)),14)-datenum('00:00:00',14);
        day=str2num(str((end-15):(end-14)));
        month=(str((end-19):(end-16)));
        switch deblank(month)
            case 'jan'
                mn=1;
            case 'feb'
                mn=2;
            case 'mar'
                mn=3;
            case 'apr'
                mn=4;
            case 'may'
                mn=5;
            case 'jun'
                mn=6;
            case 'jul'
                mn=7;
            case 'aug'
                mn=8;
            case 'sep'
                mn=9;
            case 'oct'
                mn=10;
            case 'nov'
                mn=11;
            case 'dec'
                mn=12;
        end
        
        tabs=tm+datenum(year,mn,day,0,0,0);
        
    end
end
