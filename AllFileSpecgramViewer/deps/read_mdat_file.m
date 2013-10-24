%%%%%%read_mdat_file.m%%%%%%%%%%%
% [x,fparms]=read_mdat_file(dirname,fname,tstart,nsec,calibrate_acoustic)
% Read in a DAT file.
%  WARNING!  LOG file must be in same location.
% dirname: directory where DAT and LOG file are located.
% fname, name of *DAT file.  MUST include DAT extension!
% fs, sampling frequency in Hz.
% tstart-datenumber of desired start time.  If zero, data are read from
%   start of file
% nsec-number of seconds to be pulled.
% calibrate_acoustic: if 1, then output in terms of acoustics, if 0 output in terms of voltage
%       if 2, output raw count.
%       Default is raw count
% Output of fparms:
%   tfs: datenumber of the start of the file
%   tfe: datenumber of end of file
%   tiltx, tilty:  tilt of array
%   fs:
%   nc: number of channels
%   nstart:  number of time samples into MDAT collected (
%   nsamples:   number of time samples read
% Version 1, Sept 30, 2006-read in raw data
% Version 2, Oct. 28 2006- correct for time lost from serial port
% acquisition--for this deployment every five minutes acoustic data not
% acquired for 4 seconds (plus some time to power up).  Thus I use
%   a drift formula of 48 sec/hour..


function [x,fparms]=read_mdat_file(dirname,fname,tstart,ns,calibrate_acoustic)

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
[tfs,tfe,fs,nc,tiltx,tilty,sens]=load_mdat_header(dirname,fname);

fparms=struct('tfs',tfs,'tfe',tfe,'fs',fs,'nc',nc,'tiltx',tiltx,'tilty',tilty);
fid=fopen([fname],'r','ieee-be');
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

end

function [tfs,tfe,fs,nc,tiltx,tilty,sensitivity]=load_mdat_header(dirname,fname)
cd(dirname);
fpre=findstr(lower(fname),'.mdat')-1;
log_name=[fname(1:fpre) '.LOG'];
fid=fopen(log_name,'r');
sensitivity=-160; %dB re 1uPa/V

while 1
    tline=fgetl(fid);
    if ~ischar(tline),break,end
    if findstr(tline,'Channels:')>0,
        Idot=1+findstr(tline,':');
        nc=str2num(tline(Idot:end));
    elseif findstr(tline,'Sample rate:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'Hz')-1;
        fs=str2num(tline(Idot:Iend));
    elseif findstr(tline,'Sensitivity:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'dB')-1;
        sensitivity=str2num(tline(Idot:Iend));
    elseif findstr(tline,'Start time')>0,
        Idot=findstr(tline,'=');
        tfs=parse_date(tline(Idot:end));
        tline=fgetl(fid);
        Iend=findstr(tline,'seconds')-1;
        tfs=tfs+datenum(0,0,0,0,0,str2num(tline(1:Iend)));
    elseif findstr(tline,'Stop time')>0,
        Idot=findstr(tline,'=');
        tfe=parse_date(tline(Idot:end));
        tline=fgetl(fid);
        Iend=findstr(tline,'seconds')-1;
        tfe=tfe+datenum(0,0,0,0,0,str2num(tline(1:Iend)));
    elseif findstr(tline,'Tilt-X')>0,
        Idot=findstr(tline,':')+1;
        Iend=findstr(tline,'Tilt-Y')-1;
        tiltx=str2num(tline(Idot:Iend));
        
        Idot=Iend+8;
        Iend=findstr(tline,'degrees')-1;
        tilty=str2num(tline(Idot:Iend));
        
    end
    
end

fclose(fid);
    function tabs=parse_date(str)
        
        Istart=findstr(str,'=');
        year=str2num(str((end-4):end));
        tm=datenum(str((end-13):(end-5)),14)-datenum('00:00:00',14);
        day=str2num(str((end-15):(end-14)));
        month=(str((end-19):(end-16)));
        switch deblank(month),
            case 'Jan'
                mn=1;
            case 'Feb'
                mn=2;
            case 'Mar'
                mn=3;
            case 'Apr'
                mn=4;
            case 'May'
                mn=5;
            case 'Jun'
                mn=6;
            case 'Jul'
                mn=7;
            case 'Aug'
                mn=8;
            case 'Sep'
                mn=9;
            case 'Oct'
                mn=10;
            case 'Nov'
                mn=11;
            case 'Dec'
                mn=12;
        end
        
        tabs=tm+datenum(year,mn,day,0,0,0);
        
    end
end
