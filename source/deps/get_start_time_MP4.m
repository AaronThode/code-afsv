function tmin=get_start_time_MP4(mydir,myfile,info)


tmin=datenum([1970 1 1 0 0 0]);
switch upper(myfile)
    case 'GP070503.MP4'
        tmin=datenum([2015 1 28 15 40-11 4-38]);
    case 'GP080503.MP4'
        tmin=datenum([2015 1 28 15 40 4]);
    case 'GP090503.MP4'
        tmin=datenum([2015 1 28 15 40+11 4+38]);
    case 'GP180496.MP4'
        tmin=datenum([2015 1 17 17 46 24]);
end

end
