function [x, tfs, tfe, fs]    =...
        read_dat_file(file_name, fs, tstart, ns, units_voltage, sens, raw)
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
%   raw             -   t/f, weather output should be raw data or scaled
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

%   Check that associated log file exists
log_file    =   fullfile(pathstr, [name '.log']);
if ~exist(log_file, 'file')
    log_file    =   fullfile(pathstr, [name '.LOG']); % check upper case ext
    if ~exist(log_file, 'file')
        error('%s does not exist in same directory as %s',...
            log_file,                               file_name); 
    end
end

fid	=	fopen(log_file, 'r');

line.fs     =   [];
tfs         =   0;    tfe=0;
sec_flag    =   0;

while 1
    tline	=	fgetl(fid);
    
    if ~ischar(tline), break, end
    %   ??  Date parsing can be done with inbuilt matlab functions
    if      strfind(tline, 'Start time')
        line.start  =   tline;
        tfs         =   parse_date(line.start);
        sec_flag    =   1;
    elseif  strfind(tline, 'Stop time')
        line.end    =   tline;
        tfe         =   parse_date(line.end);
        sec_flag    =   2;
    elseif  strfind(tline, 'Sample rate:')>0,
        Idot        =   1 + strfind(tline,':');
        Iend        =   strfind(tline,'Hz')-1;
        if isempty(Iend)
            Iend    =   strfind(tline,'0');
        end
        fs  =   str2double(tline(Idot:Iend(end)));
    elseif  strfind(tline,'Sensitivity:')>0,
        Idot        =   1 + strfind(tline,':');
        Iend        =   strfind(tline,'dB')-1;
        sens        =   str2double(tline(Idot:Iend));
        
    elseif  strfind(tline,'seconds')
        Is          =   strfind(tline,'seconds')-1;
        secs        =   str2double(tline(1:Is));
        if sec_flag == 1
            tfs     =   tfs + datenum(0,0,0,0,0,secs);
        elseif sec_flag == 2
            tfe     =   tfe + datenum(0,0,0,0,0,secs);
        end
        sec_flag = 0;
    end
end
fclose(fid);

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
fid     =   fopen(file_name,'r','ieee-be');
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



%%  ??  This can be done directly with datenum
function tabs=parse_date(str)

%Istart=strfind(str,'=');
year=str2double(str((end-4):end));
tm=datenum(str((end-13):(end-5)),14)-datenum('00:00:00',14);
day=str2double(str((end-15):(end-14)));
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
