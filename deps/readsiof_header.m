% rdsiof
% matlab script to read and display sample of sio file.
% whose name is in fn
% head=readsiof_header(fn);
function head=readsiof_header(fn);

%fid = fopen(char(fn),'r','ieee-le');
fid = fopen(char(fn),'r','ieee-be'); %for s1sT
%fid=fopen(char(fn),'r');
head.dlabel(1:6) = fread(fid,6,'uchar');
head.ts = fread(fid,1,'uchar');
fillx(1:57) = fread(fid,57,'uint8');
%fseek(fid,59,0);
%[ctbcx ctecx tdriftx tsampx] = fread(fid,4,'double')
a = fread(fid,4,'double');
head.ctbc=a(1);
head.ctec=a(2);
head.tdrift=a(3);
head.Fs=a(4);
fclose(fid);