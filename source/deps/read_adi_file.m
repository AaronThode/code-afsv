function [x,tfs,tfe,fs,x_raw]=read_adi_file(dirname,fname,fs,tstart,ns,units_voltage,sens)
%%%%%%read_adi_file.m%%%%%%%%%%%
% [x,tfs,tfe,fs,x_raw]=read_adi_file(dirname,fname,fs,tstart,nsec,units_voltage)
% Read in a ADIOS file, which is little endian instead of big endian, and different calibration.
%  WARNING!  LOG file must be in same location.
% dirname: directory where DAT and LOG file are located.
% fname, name of *ADI file.  MUST include ADI extension!
% fs, sampling frequency in Hz.
% tstart-datenumber of desired start time.  If zero, data are read from
%   start of file
% nsec-number of seconds to be pulled.
% units_voltage: if 1, then output in terms of voltage, otherwise output in terms of uPa.
%       Default acoustic uPa
% Output:
%   tfs: datenumber of the start of the file
%   tfe: datenumber of the end of the file
% Created Feb. 17, 2012

if ~exist('units_voltage'),
    units_voltage=0;
end
%fs=50000;
x_raw=[];
VREF=3.3;
nc=1;
if ~exist('sens')
    sens=-172; %dB re 1uPa/V
end
%Get start time of file from .log file...
mydir=pwd;
cd(dirname);
fpre=findstr(fname,'.ADI')-1;
log_name=[fname(1:fpre) '.log'];
if ~exist(log_name)
    log_name=[fname(1:fpre) '.LOG'];
     if ~exist(log_name)
        error('read_adi_file: %s does not exist in same directory as %s',log_name,fname); 
     end
end
fid=fopen(log_name,'r');

line.fs=[];
tfs=0;tfe=0;
sec_flag=0;
sensitivity_flag=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if findstr(tline,'Start time')
        line.start=tline;
        tfs=parse_date(line.start);
        sec_flag=1;
    elseif findstr(tline,'End time')
        line.end=tline;
        tfe=parse_date(line.end);
        sec_flag=2;
    elseif findstr(tline,'Sample rate:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'Hz')-1;
        if isempty(Iend)
            Iend=findstr(tline,'0');
        end
        fs=str2num(tline(Idot:Iend(end)));
        
    elseif findstr(tline,'Sensitivity:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'dB')-1;
        sens=str2num(tline(Idot:Iend));
        sensitivity_flag=1;
    elseif findstr(tline,'seconds')
        Is=findstr(tline,'seconds')-1;
        secs=str2num(tline(1:Is));
        if sec_flag==1
            tfs=tfs+datenum(0,0,0,0,0,secs);
        elseif sec_flag==2
            tfe=tfe+datenum(0,0,0,0,0,secs);
        end
        sec_flag=0;
    end
end
fclose(fid);
np=ns*fs;


if tfs==0|tfe==0
    disp(sprintf('Could not understand %s',log_name));
end
if sensitivity_flag==0
    disp(sprintf('LOG file did not have sensitivity; using a default value of %6.2f',sens));
end

try
    if tstart>tfs,
        offset=datevec(tstart-tfs);
        offset=(offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file
        %disp(sprintf('%i second offset from file start',offset));
    elseif tstart>0
        %error('cannot read this time from this file');
        disp('desired start is before file start, setting offset to 0');
        offset=0;
    else
        disp('Tstart less than or equal to zero, will interpret tstart as elasped time in seconds');
        offset=abs(tstart);
    end
catch
    disp('File name does not contain date, will interpret tstart as elasped time in seconds');
    %keyboard;
    offset=tstart;
end

scale=(1/20)*VREF/(2^16-1);  % Tested with input signal 50 mv pk-pk at 1 kHz (50MV1000.ADI)
fid=fopen(fname,'r','ieee-le');

cd(mydir);
fskip=2*nc*round(fs*offset);
res=fseek(fid,fskip,-1);
if res<0,
    keyboard;
end
%disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));
x=fread(fid,[nc np],'uint16');
fclose(fid);
if nargout>4
    x_raw=x;
end
x=scale*x;


%Crude acoustic level calibration
if units_voltage~=1,
    disp('Acoustic calibration executed');
    x=x.*(10^(-sens/20));
end
end


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

