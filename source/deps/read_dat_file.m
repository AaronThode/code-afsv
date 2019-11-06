function [x, tfs, tfe, fs,fparms]    =...
    read_dat_file(file_name, tstart, ns, units_voltage, sens, raw)
%READ_DAT_FILE     Read in a DAT file
%
%   WARNING!  LOG file must be in same location.
%
% Function form:
%  function [x, tfs, tfe, fs]    =...
%        read_dat_file(file_name, fs, tstart, ns, units_voltage, sens, raw)
%
% Input Parameters:
%   fname           -   name of *DAT file.  MUST include DAT extension and full pathname!
%   raw             -   t/f, whether output should be raw data or scaled
%   fs              -   sampling frequency in Hz.
%   tstart          -   datenumber of desired start time.
%                           If zero, data are read from start of file
%   nsec            -   number of seconds to be pulled.
%   units_voltage   -   if 1, then output in terms of voltage,otherwise
%                            output in terms of uPa. Default acoustic uPa
%   sens            -   Relative scaling factor for acoustic units
%   raw             -   If 1, output in terms of raw native units of file.
%
% Output Parameters:
%   x               -   Scaled data vector
%   tfs             -   datenumber of the start of the file
%   tfe             -   datenumber of the end of the file
%   fs              -   Sampling frequency (in case overridden internally)
%
%
% Revision History:
%   Version 1, Sept 30, 2006-read in raw data
%   Version 2, Oct. 28 2006- correct for time lost from serial port
%       acquisition--for this deployment every five minutes acoustic data
%       not acquired for 4 seconds (plus some time to power up).
%       Thus I use a drift formula of 48 sec/hour..
%
%   Modified by Jit Sarkar
%   2012/12/04
%       Header comments adjusted/formatted for readability
%   2012/12/17
%       Modified parameter checking, adding ability to read whole file
%       Modified all calls to "exist" function to be explicit, e.g.
%       searching for variables, or files
%       Changing directories just to find log file unecessary
%       Applied code-analyzer suggested fixes
%       Adjusted code formatting/spacing for readability




%%  Input parameter checking, and default values
%   dB re 1uPa/V
if ~exist('sens', 'var') || isempty(sens),      sens    =	-172;	end
if ~exist('units_voltage', 'var') || isempty(units_voltage),
    units_voltage   =   0;      end
if ~exist('fs', 'var'),                         fs      =   [];     end
if ~exist('ns', 'var'),                         ns      =   [];     end
if ~exist('tstart', 'var') || isempty(tstart),  tstart  =   0;      end
if ~exist('raw', 'var') || isempty(raw),        raw     =   false;  end


%%  Constant values and preassigned variables

VREF    =   2.5;
nc      =   1;

%   An input voltage of 0.05 V p-p produces a 1.2 V p-p at A/D
scale   =   (0.038/1.0)*VREF/(2^16-1);


%%  Parse input file name and associated log file
[pathstr, name, ~]    =   fileparts(file_name);
 pathstr=strtrim(pathstr);
 name=strtrim(name);
%   Check that associated log file exists

log_file    =   fullfile(pathstr, [name '.log']);
%Create full path name

if ~exist(log_file, 'file')
    log_file    =   fullfile(pathstr, [name '.LOG']); % check upper case ext
    if ~exist(log_file, 'file')
        error('%s does not exist in same directory as %s',...
            log_file, file_name);
    end
end

%[tfs,tfe,fs,nc,sens,Igood,geom,synch]=load_mdat_header(fname);
[tfs,tfe,fs,nc,sens0,Igood,geom,synch]=load_dat_header(log_file);

if ~isempty(sens0)
    sens=sens0;
end
fparms=struct('tfs',tfs,'tfe',tfe,'fs',fs,'nc',nc,'geom',geom,'synch',synch,'sensitivity',sens,'Igood',Igood);


if (tfs==0 || tfe==0)
    fprintf('Could not understand %s',log_file);
end


%%  Parse date/time values for reading binary file
if  isempty(fs)
    error('Sampling rate not provided, and not found in log file');
end
if  isempty(ns)     %   assume whole file is requested
    np  =   inf;
else
    np  =   ns * fs;
end


try
    if tstart > tfs
        offset  =   etime(datevec(tstart), datevec(tfs));
        %        offset  =   (offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file
        %disp(sprintf('%i second offset from file start',offset));
    elseif tstart > 0
        %error('cannot read this time from this file');
        disp('desired start is before file start, setting offset to 0');
        offset  =   0;
    else
        disp('Tstart less than or equal to zero, will interpret tstart as elasped time in seconds');
        offset  =   abs(tstart);
    end
catch    %#ok<CTCH>
    disp('File name does not contain date, will interpret tstart as elasped time in seconds');
    %keyboard;
    offset  =   tstart;
end

%%  Open binary data file and read in requeste points
fid     =   fopen(strtrim(file_name),'r','ieee-be');
fskip   =   2*nc*round(fs*offset);
res     =   fseek(fid,fskip,-1);
if res<0,
    keyboard;
end
%disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));

%   Default to raw output
x   =	fread(fid,[nc np],'*uint16');
fclose(fid);

if ~raw
    x   =   scale*double(x);
    
    %Crude acoustic level calibration
    if units_voltage~=1,
        disp('Acoustic calibration executed');
        x   =   x.*(10^(-sens/20));
    end
else
    double(x);
end

end

function [tfs,tfe,fs,nc,sens,Igood,geom,synch]=load_dat_header(log_file)
fid	=	fopen(log_file, 'r');

line.fs     =   [];
tfs         =   0;    tfe=0;
sec_flag    =   0;
nc=1;  %Channels in DAT file are always 1!
Igood=[];
synch=struct('file',[],'time',[],'offset',[],'drift',[]);
geom=struct('rd',[],'D',[],'spacing',[],'tiltx',[],'tilty',[]);
sens=[];

while 1
    tline	=	fgetl(fid);
    if ~ischar(tline),break,end
    tline_org=tline;
    tline=lower(tline);
    if ~ischar(tline), break, end
    %   ??  Date parsing can be done with inbuilt matlab functions
    if      strfind(tline, 'start time')
        line.start  =   tline;
        tfs         =   parse_date(line.start);
        sec_flag    =   1;
    elseif  strfind(tline, 'stop time')
        line.end    =   tline;
        tfe         =   parse_date(line.end);
        sec_flag    =   2;
    elseif strfind(tline,'synch time')>0
        line.synch=tline;
        synch.time=parse_date(line.synch);
        sec_flag=3;
    elseif  strfind(tline, 'sample rate:')>0,
        Idot        =   1 + strfind(tline,':');
        Iend        =   strfind(tline,'hz')-1;
        if isempty(Iend)
            Iend    =   strfind(tline,'0');
        end
        fs  =   str2double(tline(Idot:Iend(end)));
        
    elseif  strfind(tline,'sensitivity:')>0,
        Idot        =   1 + strfind(tline,':');
        Iend        =   strfind(tline,'db')-1;
        sens        =   str2double(tline(Idot:Iend));
    elseif strfind(tline,'time offset:')>0
        Idot=1+strfind(tline,':');
        Iend=strfind(tline,'sec')-1;
        synch.offset=eval(tline(Idot:Iend));
    elseif strfind(tline,'time drift:')>0
        Idot=1+strfind(tline,':');
        Iend=strfind(tline,'ms')-1;
        synch.drift=eval(tline(Idot:Iend));
    elseif  strfind(tline,'seconds')
        Is          =   strfind(tline,'seconds')-1;
        secs        =   str2double(tline(1:Is));
        if sec_flag == 1
            tfs     =   tfs + datenum(0,0,0,0,0,secs);
        elseif sec_flag == 2
            tfe     =   tfe + datenum(0,0,0,0,0,secs);
        elseif sec_flag ==3
            synch.time     =   synch.time + datenum(0,0,0,0,0,secs);
        end
        sec_flag = 0;
    elseif strfind(tline,'channels:')>0
        Idot=1+strfind(tline,':');
        nc=str2num(tline(Idot:end));
        
    elseif strfind(tline,'tilt-x')>0,
        Idot=strfind(tline,':')+1;
        Iend=strfind(tline,'tilt-y')-1;
        geom.tiltx=str2num(tline(Idot:Iend));
        Idot=Iend+8;
        Iend=strfind(tline,'degrees')-1;
        geom.tilty=str2num(tline(Idot:Iend));
    elseif strfind(tline,'synch file:')>0
        Idot=1+strfind(tline,':');
        synch.file=(tline_org(Idot:end));
        
    elseif strfind(tline,'channel depths')>0
        Idot=1+strfind(tline,':');
        geom.rd=str2num(tline(Idot:end));
    elseif strfind(tline,'channel spacing')>0
        Idot=1+strfind(tline,':');
        geom.spacing=str2num(tline(Idot:end));
    elseif strfind(tline,'water depth:')>0
        Idot=1+strfind(tline,':');
        Iend=strfind(tline,'m')-1;
        geom.D=str2double(tline(Idot:Iend));
    elseif strfind(tline,'channel quality:')>0
        Idot=1+strfind(tline,':');
        Igood=str2num(tline(Idot:end));
    end
    
end
fclose(fid);

end %header


%%  ??  This can be done directly with datenum
function tabs=parse_date(str)

%Istart=strfind(str,'=');
year=str2double(str((end-4):end));
tm=datenum(str((end-13):(end-5)),14)-datenum('00:00:00',14);
day=str2double(str((end-15):(end-14)));
month=(str((end-19):(end-16)));
switch deblank(month),
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
