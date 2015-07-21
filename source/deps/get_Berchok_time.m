function [SUDAR_true,tmin,tmax]=get_Berchok_time(mydir,myfile,Nsamples,Fs) %Check whether a sUDAR file exists
SUDAR_true=false;
tmin=[];
tmax=[];

[pathstr,fname,extt] = fileparts(myfile);
fname_log=fullfile(mydir,[fname  extt]);
fname_found=dir(fname_log);

if isempty(fname_found)
    tmin=[];
    tmax=[];
    
    return
end
%AU-BS02b-121014-192730
Idot=max(strfind(myfile,'.'));
Idash=strfind(myfile,'-');
start_str=myfile((Idash(end-1)+1):(Idash(end)-1));
end_str=myfile((Idash(end)+1):(Idot-1));

try
    %start_str1=eval(start_str);
    %end_str1=eval(end_str);
    yr=eval(start_str(1:2));
    mth=eval(start_str(3:4));
    day=eval(start_str(5:6));
    hr=eval(end_str(1:2));
    minn=eval(end_str(3:4));
    sec=eval(end_str(5:6));
    
    tmin=datenum(2000+yr,mth,day,hr,minn,sec);
    tmax=tmin+datenum(0,0,0,0,0,Nsamples/Fs);  %Five min long?
catch
    return
end


SUDAR_true=true;

end

function tm=parse_data(tline)
left=min(strfind(tline,'>'))+1;

right=max(strfind(tline,'</'))-1;
tm=(tline(left:right));

end

