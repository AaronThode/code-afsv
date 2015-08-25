function [x, tfs, tfe, fs,fparms]    =...
    read_synchronized_dat_file(file_name, tstart, ns, units_voltage, sens, raw)
%READ_DAT_FILE     Read in a DAT file
%
%   WARNING!  LOG file must be in same location.
%
% Function form:
%  function [x, tfs, tfe, fs]    =...
%        read_synchronized_dat_file(file_name, fs, tstart, ns, units_voltage, sens, raw)
%
% Input Parameters:
%   fname           -   name of *DAT file.  MUST include DAT extension and full pathname!
%   raw             -   t/f, weather output should be raw data or scaled
%   tstart          -   datenumber of desired start time.
%                           If zero, data are read from start of file
%   ns            -   number of seconds to be pulled.
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




%%  Input parameter checking, and default values
%   dB re 1uPa/V
if ~exist('sens', 'var') || isempty(sens),      sens    =	-172;	end
if ~exist('units_voltage', 'var') || isempty(units_voltage),
    units_voltage   =   0;      end
if ~exist('ns', 'var'),                         ns      =   [];     end
if ~exist('tstart', 'var') || isempty(tstart),  tstart  =   0;      end
if ~exist('raw', 'var') || isempty(raw),        raw     =   false;  end


%Try to get sample frequency from from file...
[~,~,~,fs]=read_dat_file(file_name,-1,1,0);

if isempty(fs)
    fs=input('Enter sampling rate in Hz:');
end

%Set up directory navigation...
[dirname,fname,extt] = fileparts(file_name);
fname=[fname extt];
dirname=dirname(~isspace(dirname));
thisdir=pwd;
cd(dirname)

%Read in basic file
[xbot,tfs,tfe,fs,fparms_bot]=read_dat_file(fname,tstart,ns,units_voltage); %output in uPa
fparms=fparms_bot;

if isempty(fparms.synch.file)
    x=xbot;
    if size(x,2)>1
        x=x';
    end
    cd(thisdir);
    return
end

%[xbot,fparms_bot]=read_mdat_file(fname,tdate_start,tlen,1);


if tstart<0
    tstart=fparms_bot.synch.time;
end

elapsed_time=tstart-fparms_bot.synch.time;
if elapsed_time>0
    etime=datevec(elapsed_time);
    nsec=3600*etime(4)+60*etime(5)+etime(6);
else
    etime=datevec(-elapsed_time);
    nsec=3600*etime(4)+60*etime(5)+etime(6);
    nsec=-nsec;
end

toffset=fparms_bot.synch.offset+fparms_bot.synch.drift*nsec/(3600*1000);
fprintf('Toffset is %5.6f msec\n',toffset*1000);
try  %Look for corresponding file
    [xtop,~,~,~,fparms_top]=read_dat_file(fparms.synch.file,tstart+datenum(0,0,0,0,0,toffset),ns,units_voltage);
    
    x=[xbot' xtop'];  %Start with bottom element, move to top.
    
    rd=[fparms_bot.geom.rd fparms_top.geom.rd];
    %spacing=[fparms_bot.geom.spacing fparms_top.geom.spacing];
    Igood=[fparms_bot.Igood fparms_top.Igood];
catch
    Igood=fparms_bot.Igood;
    x=xbot;
    rd=fparms_bot.geom.rd;
end
%Include only valid channels..
rd=rd(Igood>0);
Ichan=1:size(x,2);
x=x(:,Igood>0);
%Ichan=Ichan(Igood>0);

[rd,Iorder]=sort(rd);  %Shallowest element is first...
%rd=fliplr(rd);
%Igood=Igood(Iorder);
x=x(:,Iorder);
%Igood=find(Igood(Iorder)>0);
spacing=[0 diff(rd)];

fparms=fparms_bot;
fparms.geom.rd=rd;
fparms.geom.Igood=Igood;
fparms.geom.Iorder=Iorder;
fparms.geom.spacing=spacing;

%x=flipud(x);  %Shallowest element now first element, matching head.rd
cd(thisdir);

%%Short-circuit data load with synthesized file...

end
