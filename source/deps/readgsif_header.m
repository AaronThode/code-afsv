% head=readgsif_header(fn);
%  Inputs header information from a GSI file
%  fn: file name string (includes .gsi)
%  head: structure containing:
%		dlabel
%		nc
%		ctbc,ctec,tdrift
%		Fs, UTMX,UTMY,depth,UTMZone
%		brefa, ts,tabs_start,tabs_end
function head=readgsif_header(fn)

%fid = fopen(char(fn),'r','ieee-le');
fid = fopen(char(fn),'r','ieee-be'); %for s1sT
head.dlabel(1:10) = char(fread(fid,10,'uchar'));
%head.nc=fread(fid,1,'uchar');
head.contents=char(fread(fid,4,'uchar'));
head.nc=length(find(double(head.contents)>34));

fillx(1:50) = fread(fid,50,'uint8');

a = fread(fid,9,'double');
head.ctbc=a(1);
head.ctec=a(2);
head.tdrift=a(3);
head.Fs=a(4);
head.UTMX=a(5);
head.UTMY=a(6);
head.depth=a(7);
head.UTMZone=a(8);
head.brefa=a(9);
head.ts = fread(fid,1,'uchar');
head.tabs_start=datenum(1970,1,1,0,0,head.ctbc);
head.tabs_end=datenum(1970,1,1,0,0,head.ctec);

fclose(fid);
