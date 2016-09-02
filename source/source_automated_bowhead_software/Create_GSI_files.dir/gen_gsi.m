% gen_gsi.m
%This version of Gen_gsi.m is the same as in V1, except it is
%specific for the data taken for the Shell arrays in 2007.
%These DASARs had four data channels (Omni, cosine, sine, and vertical).
%The vertical channels were never of much use.
%This version of gen_gsi will read the four-channel data,
%and convert it to the 3-channel version standard in 2008 and
%earlier years, omitting the vertical channel.
% Matlab script to generate one 3-channel time series file for
% Scripps automated call detection processing.
% It is called from Supergen_gsi, which sets inputs.
% Inputs, set before invoking gen_gsi, are:
% yearsite (set by rdsumry program)
% datapath (set by rdsumry program)
% id,DasarID   dasar number and label   (set by invoking program)
% mtbc Matlab start time vector  "
% mtec Matlab end time vector
%
% Input is DASAR raw time series file with continuous 8-byte samples
% consisting of pressure, x-vel, y-vel, z-vel
% Output file is continuous pressure, x-vel, y-vel time series covering a specified time.
% output file for 1 day is about 172.8 MB
% Output has a header containing
% 10-character DASAR ID consisting of year,site, and dasar letter (eg D07_1a)
% 1 byte specifying number of a/d channels (3)
% 4-characters specifying  P=pressure   X=x component  Y=y component blank(or Z)
% 1-char specifying use:T=use all;F=don't use;D=detece,but no loc
% 48 bytes reserved
% 64-byte floating ctbase absolute ctime of first sample
% 64-byte floating ctend  absolute ctime of last sample
% 64-byte floating tdrift time drift of DASAR clock in sec/day
% 64-byte floating nominal sample rate(not including clock drift)
% 64-byte floating UTM X in meters of DASAR Location
% 64-byte floating UTM Y
% 64-byte floating depth, meters
% 64-byte floating UTM Zone
% 64-byte floating brefa, grid bearing in degrees of DASAR cosine axis
% 376 bytes reserved
%
%-----------------------------------------------
% some terminology:
% 1 adcr is one 16-bit a/d converter reading on one channel =2bytes
% 1 sample is a group of 3 or 4 adcr s taken at one time
% chpis is number of adcrs per input sample (=4 for 2007; =3 for 2008
% chpos is number of adcrs per output sample (=3 for 2007 and 2008)
% chpis = 4;  %adcrs per input sample
% chpos = 3;  %adcrs per output sample
% rdsumry08 must have been run.
%
%-----------------------------------------------
% check that output disk has folders ready.

function tba=gen_gsi(did,utmzone,paths,mtbc,mtec,chpsin,chpsout,buffer_samples)

DasarID=char(did.lbl);
tsamp = 1000.; %Sampling rate in Hz, without clock drift
bbo =   33554432;  % this if processing from full data files
%bpsin = 8;  %Bytes per sampe of time (i.e 8 for four-channel 16 bit raw data
%bpsout = 6;

bpsin = chpsin*2;  %Bytes per sampe of time (i.e 8 for four-channel 16 bit raw data)
bpsout = chpsout*2;

%

% pre-allocate output buffer
% this is the path to the raw data files:
ctbc = mat2c_tm(mtbc);      %convert to ctime
ctec = mat2c_tm(mtec);      % convert to c-time
%tseries = 'P';

%
% compute starting input byte address in raw files
%ctsp = ctbc;
tbai=ctime2byteoffset(ctbc,did,tsamp,bpsin,bbo);
tbaf=ctime2byteoffset(ctec,did,tsamp,bpsin,bbo);

totalsamples = (tbaf-tbai)/bpsin +1;
totalbuffers = totalsamples/buffer_samples;  %in output file
wholebuffers = floor(totalbuffers);
fractbuffer = totalbuffers -wholebuffers;
fractsamps = fractbuffer*buffer_samples;  %Number of samples left in remainder
%
% we will fill a buffer_samples buffer totalbuffers+1 times and write it to the
% output file.
% form output file name and open it
%ofd = [paths.dasarfolder '/'];
%ofn = char(strcat(ofd, DasarID,'T', datestr(mtbc,30), '.gsi'));
ofn=sprintf('%s/%sT%s.gsi',paths.dasarfolder,DasarID,datestr(mtbc,30));
ofid = fopen( ofn,'w','ieee-be');
%
% write header sector
writeGSIheader(ofid,did,ctbc,ctec,utmzone,tsamp,chpsout);

tba = tbai;


for iob = 1:wholebuffers+1
    disp(sprintf('%6.1f percent done...',100*iob/(wholebuffers+1)));
    if iob==wholebuffers+1
        nw = round(fractsamps);
    else
        nw = buffer_samples;
    end
    
    tba=ReadRAWDataWriteGsi(ofid,tba,buffer_samples,chpsin,chpsout,paths,nw);
    
end
fclose(ofid);
end



