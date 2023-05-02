function [SUDAR_true,tmin,tmax,Fs,cal_dB,tmin_UTC,tmax_UTC]=get_SUDAR_time(mydir,myfile)
%function [SUDAR_true,tmin,tmax,Fs,cal_dB,tmin_UTC,tmax_UTC]=get_SUDAR_time(mydir,myfile)
%Check whether a SUDAR (SoundTrap) file exists by looking for log.xml file

SUDAR_true=false;

[pathstr,fname,~] = fileparts(myfile);
fname_log=fullfile(mydir,[fname '.log.xml']);
fname_found=dir(fname_log);

if isempty(fname_found)
    tmin=[];
    tmax=[];
    Fs=[];
    return
end
SUDAR_true=true;

fid=fopen(fname_log,'r');
while 1
    tline	=	fgetl(fid);
    
    if ~ischar(tline), break, end
    %   ??  Date parsing can be done with inbuilt matlab functions
    if contains(tline,'AUDIO Gain="Low"')
        cal_dB=185.8;
    elseif contains(tline,'AUDIO Gain="High"')
        cal_dB=173.1;
    elseif      strfind(tline, 'SamplingStartTimeLocal=')>0
        %tmin         =  datenum(parse_time(tline));  %%Works for me
        datte=parse_time(tline);
        %tmin=datenum(datte(1:10))+datenum(datte(12:end),2019)-datenum(2019,1,1,0,0,0);
        % %% Fixing a soundtrap with the wrong year.
        
        tmin=datenumm(datte);
        
        
    elseif  strfind(tline, 'SamplingStopTimeLocal=')>0
        %tmax         =   datenum(parse_time(tline));
        datte=parse_time(tline);
        %tmax=datenum(datte(1:10))+datenum(datte(12:end),2019)-datenum(2019,1,1,0,0,0);
        tmax=datenumm(datte);
    elseif  strfind(tline, 'SamplingStartTimeUTC=')>0
        %tmax         =   datenum(parse_time(tline));
        datte=parse_time(tline);
        %tmax=datenum(datte(1:10))+datenum(datte(12:end),2019)-datenum(2019,1,1,0,0,0);
        tmin_UTC=datenumm(datte);
        
    elseif  strfind(tline, 'SamplingStopTimeUTC=')>0
        %tmax         =   datenum(parse_time(tline));
        datte=parse_time(tline);
        %tmax=datenum(datte(1:10))+datenum(datte(12:end),2019)-datenum(2019,1,1,0,0,0);
        tmax_UTC=datenumm(datte);
        
    elseif  strfind(tline, '<FS>')>0
        Fs  =   str2double(parse_data(tline));
        
        %Fix a bug
        %Fs=Fs/2;
    end
    %     elseif  strfind(tline,'Sensitivity:')>0,
    %         Idot        =   1 + strfind(tline,':');
    %         Iend        =   strfind(tline,'dB')-1;
    %         sens        =   str2double(tline(Idot:Iend));
    %
    
end
fclose(fid);

%Option to convert to local time if not accurate..
%Convert to local time using filename as a clue...
%Idot=findstr(fname,'.')+7;
%hr=str2num(fname(Idot:(Idot+1)));

%temp=datevec(tmin);
%time_zone=temp(4)-hr;  %To subtract from all times
%temp(4)=(hr);
%tmin=datenum(temp);

%temp=datevec(tmax);
%temp(4)=temp(4)-time_zone;
%tmax=datenum(temp);

end

function tm=datenumm(datte)
%%%%Try to handle odd formats of XML datestru
try
    tm=datenum(datte);
catch
    tm=datenum(datte,'yyyy-mm-ddTHH:MM:SS');
end
end

function tm=parse_data(tline)
left=min(strfind(tline,'>'))+1;

right=max(strfind(tline,'</'))-1;
tm=(tline(left:right));

end

function tm=parse_time(tline)
left=min(strfind(tline,'"'))+1;

right=max(strfind(tline,'"'))-1;
tm=(tline(left:right));

end




