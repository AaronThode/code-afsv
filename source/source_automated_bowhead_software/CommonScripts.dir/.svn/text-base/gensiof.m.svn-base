% gensiofx8.m
% Matlab script to generate pressure time series file for Scripps
% Processing.
% Input is DASAR raw time series file with continuous 8-byte samples
% consisting of pressure, x-vel, y-vel, z-vel
% Output file is continuous pressure time series.
% output file for 1 day is about 172.8 MB
% Output has a header containing
% 6-character DASAR ID consisting of year,site, and dasar letter (eg D07_1a)
% 1-character specifying  P=pressure   X=x component  Y=y component Z=z component
% 57 bytes reserved
% 64-byte floating ctbase absolute ctime of first sample
% 64-byte floating ctend  absolute ctime of last sample
% 64-byte floating tdrift time drift of DASAR clock in sec/day
% 64-byte floating nominal sample rate(not including clock drift)
% 192 bytes reserved
% Path to raw data files is specified below
% as are DASAR ID and start and end times.
% rdsumry must have been run and set to correct site.
% This variation has same logic as gensiof, but 8 times the buffer size
% trying to speed it up
%
%-----------------------------------------------
tic   %start elapsed timer
bsize = 8192;   %samples
% this is the path to the raw data files:
prdf = '\\fin\dasar\1\shell06\raw06\b1\';
mtbc = [2006 9 21 0 0 0];  % begin copy time
ctbc = mat2c_tm(mtbc);      %convert to ctime
mtec = [2006 9 22 0 0 0];   % end copy time
ctec = mat2c_tm(mtec);      % convert to c-time
cdid = 'B1';                 % DASAR to be copied
ofd = 'c:\shell07\site07_1\SIOinput07_1\';  % output file directory
tseries = 'P';
tsamp = 1000.;
bbo =   33554432;  % this if processing from full data files
%
% get DASAR id by comparing did to ID in dsumry
for id = 1:idmax
    if strcmp(cdid,did(id).lbl)
        break
    end
end
%    
% compute starting byte address in raw files
ctsp = ctbc;
if ctdstart(id)<ctsp & ctsp<ctdend(id) & strcmp(use(id),'T')
    % compute time correction to sample clock
    dtoff=tint(id) + (ctsp-ctref(id))*tdrift(id)/86400;
    ctsamp=ctsp+dtoff;
    %compute offset from file start rounded to nearest ms
    toff=ctsamp-ctstart(id);    %seconds and fractions
    tbai=round(1000*toff)*8 +bbo   %even multiple of 8 bytes
    
end
% compute ending byte address
ctsp = ctec;
if ctdstart(id)<ctsp & ctsp<ctdend(id) & strcmp(use(id),'T')
    % compute time correction to sample clock
    dtoff=tint(id) + (ctsp-ctref(id))*tdrift(id)/86400;
    ctsamp=ctsp+dtoff;
    %compute offset from file start rounded to nearest ms
    toff=ctsamp-ctstart(id);    %seconds and fractions
    tbaf=round(1000*toff)*8 +bbo -8;   %even multiple of 8 bytes
end
totalsamples = (tbaf-tbai)/8 +1;
totalbuffers = totalsamples/bsize;  %in output file
wholebuffers = floor(totalbuffers);
fractbuffer = totalbuffers -wholebuffers;
fractsamps = fractbuffer*bsize;
%
% we will fill a bsize buffer totalbuffers+1 times and write it to the
% output file.
% form output file name and open it
ofn = [ofd 'D06_1' cdid '_' datestr(mtbc,30) '.sio'];
ofid = fopen( ofn,'w','ieee-be');
%
% write header sector
count1 = fwrite(ofid,['D06_1 ' tseries ]);
fill(1:57)=0;
prec = char('double');
count2 = fwrite(ofid,fill(1:57));
count3 = fwrite(ofid,[ctbc ctec tdrift(id) tsamp], 'double');
fill(1:416)=0;
count4 = fwrite(ofid,fill(1:416));
%
tba = tbai;
for iob = 1:totalbuffers+1
    % see if time for operator message
    % compute file and byte address of nba
    iraw=floor(tba/512000000);
    baddr=tba-512000000*iraw;
    fstr=int2str(200+iraw);
    DasarID=did(id).lbl;
    fn = char(strcat(prdf, DasarID,fstr(2:3),'.raw'));
%    disp(['baddr=' sprintf('%12.0d',baddr)])
    % number of samples to read
    K=bsize*8;
    %read time series data into array
    [fid,message] = fopen(fn,'r','ieee-be') ;
    status = fseek(fid,baddr,-1);
    dataok(id)=status~=-1;
    if(status~=0);
        ferror(fid)
    end
%if file opened and seek was ok read data.
    if status ~= -1
        % data is big-endian (Motorola)
        [x,count] = fread(fid,[4,K],'ushort'); % 8 = 4 ch x 2 bytes/sample
        %count
        fclose(fid);
         %if record is short, read rest from start of next raw file
        if count < 4*K
            iraw=iraw+1;
            fstr=int2str(200+iraw);
            fn = char(strcat(prdf,DasarID,fstr(2:3),'.raw'));
            [fid,message] = fopen(fn,'r','ieee-be');
            status=fseek(fid,0,-1);
            [x(:,length(x)+1:K),count2] = fread(fid,[4,K-length(x)],'ushort');
            fclose(fid);
        end
%        if iob == 1
%            x(1:4,1)
%            sprintf('%6.0X',x(1:4,1)),
%        end
%*************************
        if iob==1
            x(:,1:2)
        end
%*************************
        x = x';
        %New scaling is:
        sclf=2.5/65536;
        om = x(:,1);
%        cosch(:,id) = x(:,3)*sclf;
%        sinch(:,id) = x(:,2)*sclf;
%        zch(:,id) = x(:,4)*sclf;

        % strip off DC
%        om = om - round(mean(om));
%        cosch(:,id) = cosch(:,id) - mean(cosch(:,id));
%        sinch(:,id) = sinch(:,id) - mean(sinch(:,id));
%        zch(:,id) = zch(:,id) -mean(zch(:,id));
    end
% om is now just the 1-d pressure time series. Write to output
% If this is last sector, clear unused part of om buffer
    if iob==totalbuffers+1
        nw = fractsamp;
    else
        nw = bsize;
    end
    count5 = fwrite(ofid,om(1:nw),'uint16');
    tba=tba+8*bsize;
end
fclose(ofid);
toc     %stop elapsed timer
%disp(['Time to copy 1 day=' num2str(toc)])