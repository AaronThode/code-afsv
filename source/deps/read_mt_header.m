%%%%%%%%%%%read_mt_header.m%%%%%%%
% Read in mt header information
%function head=read_mt_header(fname),
function head=read_mt_header(fname)

if isempty(strfind(lower(fname),'.mt'))
    fname=[fname '.mt'];    
end
%fname='aa_Sound_2003_08120945.mt';
fid=fopen(fname,'r');
head.magicstring=read_char(fid,8);
head.totalhdrs=read_char(fid,3);
head.abbrev=read_char(fid,8);
head.stationcode=read_char(fid,3);
head.title=read_char(fid,82);
head.month=read_num(fid,3);
head.day=read_num(fid,3);
head.year=read_num(fid,5);
head.hours=read_num(fid,3);
head.minutes=read_num(fid,3);
head.seconds=read_num(fid,3);
head.msec=read_num(fid,4);
head.sampling_period=read_char(fid,15);
head.Fs=1./str2num(head.sampling_period);
head.samplebits=read_num(fid,3);
head.wordsize=read_num(fid,2);
head.typemark=read_char(fid,1);
head.swapping=read_char(fid,1);
head.signing=read_char(fid,1);
head.caltype=read_char(fid,1);
head.calmin=read_num(fid,15);
head.calmax=read_num(fid,15);
head.calunits=read_char(fid,40);
head.recordsize=read_num(fid,6);
head.sourcevers=read_char(fid,9);
head.sourcesn=read_char(fid,16);
head.tstart=datenum(head.year,head.month,head.day,head.hours, ...
    head.minutes,head.seconds+head.msec/1000);

fseek(fid,0,1);
byte_length=ftell(fid);
sec_length=byte_length/(head.wordsize*head.Fs);
head.tend=datenum(head.year,head.month,head.day,head.hours, ...
    head.minutes,head.seconds+head.msec/1000+sec_length);

fclose(fid);

function nhout=read_num(fid,N)
nhout=str2num(char(fread(fid,N,'char')'));

function chout=read_char(fid,N)
chout=char(fread(fid,N,'char')');