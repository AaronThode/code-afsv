
function tba_out=ReadRAWDataWriteGsi(ofid,tba,nsample,chpsin,chpsout,paths,nw)
%tba-byte offset desired from start of raw data files
%nsample: number of samples to read in

bpsin = chpsin*2;  %Bytes per sampe of time (i.e 8 for four-channel 16 bit raw data)
if chpsin~=chpsout
    error('ReadRAWDataWriteGsi: chpsin must equal chpsout');
    return
end
obuf = zeros(chpsout,nsample,'uint16'); %Total output buffer size in adcrs
%

iraw=floor(tba/512000000);
baddr=tba-512000000*iraw;
fstr=int2str(1000+iraw);
%fn = char(strcat(datapath,'Raw',earsite,'\',DasarID,'\', DasarID,fstr(3:4),'.raw'));
fn =sprintf('%s/%s_%s.raw',paths.raw,paths.DasarID,fstr(2:4));
%    disp(['baddr=' sprintf('%12.0d',baddr)])

% number of samples to read
%K = nsample;  %samples
%read time series data into array
[fid,message] = fopen(fn,'r','ieee-be') ;
status = fseek(fid,baddr,-1);
%dataok(id)=status~=-1;
if(status~=0);
    ferror(fid)
end
%if file opened and seek was ok read data.
if status ~= -1
    % data is big-endian (Motorola)
    [x,count] = fread(fid,[chpsin,nsample],'uint16'); % 8 = 4 ch x 2 bytes/sample
    %count
    fclose(fid);
    
    obuf(:,1:size(x,2))=x;
    
    %if record is short, read rest from start of next raw file
    if count < chpsin*nsample

        iraw=iraw+1;
        fstr=int2str(1000+iraw);
        fn =sprintf('%s/%s_%s.raw',paths.raw,paths.DasarID,fstr(2:4));
        [fid,message] = fopen(fn,'r','ieee-be');
        status=fseek(fid,0,-1);

	%AARON correction:  Since a 6-byte sample (2-btye three channels) is not
% evenly divided into 512MB, the last column of x may be unfilled 
% i.e. the first sample of the following raw file may not be the pressure sample.
	Izero=find(x(:,end)==0);
	if ~isempty(Izero)
		obuf(Izero,size(x,2))=fread(fid,[length(Izero),1],'uint16');
	end
        [obuf(:,size(x,2)+1:nsample),count2] = fread(fid,[chpsin,nsample-size(x,2)],'uint16');
        %[x(:,length(x)+1:K),count2] = fread(fid,[chpsin,K-length(x)],'uint16');
        fclose(fid);
        
    end
    
    
end

% obuf is now the 3-channel time series. Write to output
% If this is last sector, clear unused part of om buffer

count5 = fwrite(ofid,obuf(:,1:nw),'uint16');
tba_out=tba+bpsin*nsample;
end
