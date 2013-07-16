% readEnergy_header(fname)
% matlab script to read header of JAVA energy detection file
% whose name is  fname
% 
% Updated Jan 7, 2010 to incorporate channel size and a potential brefa
%   (reference bearing provided by detections derived from GSI files)

function head=readEnergy_header(fn)

%fid = fopen(char(fn),'r','ieee-le');
fid = fopen(char(fn),'r','ieee-be'); %for s1sT
%fid=fopen(char(fn),'r');
head.dlabel = char(fread(fid,38,'uchar')');
head.dlabel=head.dlabel(2:2:end);
head.Fs = fread(fid,1,'double');
head.Ndetectors = fread(fid,1,'int32');
head.Nfft = fread(fid,1,'int32');
head.dn = fread(fid,1,'int32');
head.threshold=fread(fid,1,'double');
head.MinimumDetectionTime=fread(fid,1,'float32');
head.eq_time=fread(fid,1,'float32');
for I=1:head.Ndetectors,
    head.flow(I)=fread(fid,1,'int32');
    head.fhigh(I)=fread(fid,1,'int32');
end

head.nchan=fread(fid,1,'int32');

%%Check to see if bearing information exists (i.e. were these detections
%%derived from a GSI file?

test=fread(fid,1,'int32');
if test==1
    head.brefa=fread(fid,1,'double');
else
    %keyboard
end
head.wordtype=char(fread(fid,12,'uchar')');
head.wordtype=head.wordtype(2:2:end);
fclose(fid);
