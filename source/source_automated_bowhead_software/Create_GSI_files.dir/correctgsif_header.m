% correctgsif_header
% Correct header of gsi file: change brefa, tdrift, or ctbc 
% whose name is in fn
%function correctgsif_header(fn,brefa,tdrift,ctbc)
function correctgsif_header(fn,brefa,tdrift,ctbc)

%fid = fopen(char(fn),'r','ieee-le');
fid = fopen(char(fn),'r+','ieee-be'); %for s1sT
head.dlabel(1:10) = char(fread(fid,10,'uchar'));
%head.nc=fread(fid,1,'uchar');
head.contents=char(fread(fid,4,'uchar'));
head.nc=length(deblank(head.contents'));

fillx(1:50) = fread(fid,50,'uint8');

a = fread(fid,1,'double');
head.ctbc=a(1);
if exist('ctbc')
    fseek(fid,-8,0); %Move one double back
    fwrite(fid,ctbc,'double');

end
a = fread(fid,2,'double');

head.ctec=a(1);
head.tdrift=a(2);

if exist('tdrift')
    fseek(fid,-8,0); %Move one double back
    fwrite(fid,tdrift,'double');
end
a = fread(fid,6,'double');

%fseek(fid,-8,0); %Move one double back
%brefan=fread(fid,1,'double');
if exist('brefa')
    fseek(fid,-8,0); %Move one double back
    fwrite(fid,brefa,'double');
end
%head.ts = fread(fid,1,'uchar');
%head.tabs_start=datenum(1970,1,1,0,0,head.ctbc);

fclose(fid);
