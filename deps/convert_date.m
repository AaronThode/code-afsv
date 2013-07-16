%%%%%%convert_date.m%%%%%
% Extract a datenumber from a filename string...
% function tstart=convert_date(data,delimiter),
%%get start date from file name.
%%Two types of filenames, one "MT" style the other MATLAB 'T" style
%%Start intelligent search for a timestamp in string name...


function tstart=convert_date(data,delimiter),

mt_style_flag=~isempty(findstr('Sound',data));
fname_bounds=findstr(data,delimiter);

if length(fname_bounds)==1,
    fname_bounds(2)=length(data)+1;
end
Idot=findstr(data,'.');
fname_bounds=[fname_bounds Idot(1)];

nm=[];
if mt_style_flag==1,
    Isound=findstr(data,'Sound')+6;
    yr=str2num(data(Isound:(Isound+3)));
    mo=str2num(data((Isound+5):(Isound+6)));
    dy=str2num(data((Isound+7):(Isound+8)));
    hr=str2num(data((Isound+9):(Isound+10)));
    mn=str2num(data((Isound+11):(Isound+12)));
    %sec=str2num(data((Isound+13):(Isound+14)));
    tstart=datenum(yr,mo,dy,hr,mn,0);
    
else
    try,
        
        for I=1:(length(fname_bounds)-1),
            test=data((fname_bounds(I)+1):(fname_bounds(I+1)-1));
            if length(test)==15,
                if strcmp(test(9),'T'),
                    nm=test;
                end
            end
        end
        if ~isempty(nm)
            yr=str2num(nm(1:4));
            mo=str2num(nm(5:6));
            dy=str2num(nm(7:8));
            hr=str2num(nm(10:11));
            mn=str2num(nm(12:13));
            sec=str2num(nm(14:15));
            
            tstart=datenum(yr,mo,dy,hr,mn,sec);
        end
    catch, 
        error('String does not contain valid datestr');
    end
end