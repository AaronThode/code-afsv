function str=ctime2str(ctime)

if ~isinf(ctime),
    str=datestr(datenum(1970,1,1,0,0,ctime),31);
else
    str='Inf';
end