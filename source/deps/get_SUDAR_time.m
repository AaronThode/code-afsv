function [SUDAR_true,tmin,tmax,Fs]=get_SUDAR_time(mydir,myfile) %Check whether a sUDAR file exists
  SUDAR_true=false;
          
[pathstr,fname,~] = fileparts(myfile);
fname_log=fullfile(mydir,[fname '.log.txt']);
fname_found=dir(fname_log);

if isempty(fname_found)
    tmin=[];
    tmax=[];
    return
end
SUDAR_true=true;
  
fid=fopen(fname_log,'r');
while 1
    tline	=	fgetl(fid);
    
    if ~ischar(tline), break, end
    %   ??  Date parsing can be done with inbuilt matlab functions
    if      strfind(tline, '<Sampling Start Time>')>0
        tmin         =  datenum(parse_data(tline));
        
    elseif  strfind(tline, '<Sampling Stop Time>')>0
         tmax         =   datenum(parse_data(tline));
       
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

%Convert to local time using filename as a clue...

Idot=findstr(fname,'.')+7;
hr=str2num(fname(Idot:(Idot+1)));

temp=datevec(tmin);
time_zone=temp(4)-hr;  %To subtract from all times
temp(4)=(hr);
tmin=datenum(temp);

temp=datevec(tmax);
temp(4)=temp(4)-time_zone;
tmax=datenum(temp);

end

function tm=parse_data(tline)
left=min(strfind(tline,'>'))+1;

right=max(strfind(tline,'</'))-1;
tm=(tline(left:right));

end

